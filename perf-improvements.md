# PureCLIP Performance Improvement Opportunities

Analysis of `src/` — ~7,000 lines of C++ header-only HMM pipeline.
Benchmark target: human genome (hg19, 2.9 GB), 2.9M reads per replicate, 1–2 replicates.

---

## Profile summary (estimated)

| Phase | Est. share | Parallel? | Bottleneck |
|-------|-----------|-----------|------------|
| KDE computation | **40–55 %** | ✅ OpenMP (coarse) | O(n²) brute-force convolution |
| Forward-Backward | **15–25 %** | ✅ OpenMP (per interval) | Repeated alloc/free, NaN checks |
| Emission probs | **10–15 %** | ✅ OpenMP (per interval) | `exp`/`log`/`lgamma` per position |
| GSL parameter fitting | **5–15 %** | ❌ Serial | Nelder-Mead over full dataset |
| BAM I/O + interval build | **5–10 %** | ❌ BGZF=1 thread | Serial decompression |

---

## 1. KDE computation — the dominant bottleneck

**Where:** `Observations::computeKDEs()` in `src/util.h` (lines ~397–450).

**What it does:** For every covered genomic interval (thousands per chromosome), computes a Gaussian or Epanechnikov kernel density estimate via a sliding weighted sum:

```cpp
for (unsigned t = 0; t < length(); ++t) {
    double kde = 0.0;
    for (unsigned i = max(t - w_50, 0); i < length() && i <= t + w_50; ++i) {
        kde += truncCounts[i] * kernelDensities[abs(t - i)];
    }
    kdes[t] = kde / bandwidth;
}
```

This is a **brute-force O(T × W) convolution** with window size W = 4 × bandwidth (default: 200 positions). For a contig like chr1 (249 Mb), even sparse regions trigger thousands of interval scans. There are two full passes per interval (`kdes` + `kdesN`), and this runs twice per chromosome during learning (`computeSLR` phase + `preproCoveredIntervals`) plus once during the apply phase.

**Why it's slow:**
- O(T × W) instead of O(T log T)
- Two separate KDE passes (different kernels/bandwidths) → double the work
- No SIMD vectorization (`std::abs` + float multiply in inner loop = branchy)
- Kernel weights precomputed but the code comment itself notes: *"inefficient if low genome coverage and not selected only for covered regions!!!"*

### Proposed fix: FFT-based convolution

Replace the O(T × W) windowed sum with an **FFT-based convolution** using `fftw3` or `vDSP` (Accelerate on macOS). For each interval of length T:

```
kde = IFFT(FFT(truncCounts_padded) × FFT(kernel_padded)) / bandwidth
```

This drops the inner-loop complexity from O(T × W) → O(T log T) and can leverage SIMD automatically through the FFT library. For a 249 Mb chromosome, the speedup is substantial — convolution costs go from roughly billons of operations to millions.

**Estimated impact:** 10–100× on KDE, 2–5× overall wall time reduction.

**Alternative (simpler):** If FFT integration is too invasive, add `#pragma omp simd` + restrict pointers to the inner loop and vectorize the broadcast multiply-accumulate. Even without FFT this can give 4–8× on the KDE phase via SIMD alone.

---

## 2. Forward-Backward algorithm — allocation churn + NaN guarded log-space

**Where:** `HMM::computeStatePosteriorsFB()`, `iForward()`, `iBackward()` in `src/hmm_1.h` (lines ~596–830).

**What it does:** For each covered interval, allocates `alphas_1` and `betas_1` arrays (T × 4 of `long double`), runs the log-space forward and backward passes (4 additions + `get_logSumExp_states()` per t, per k), then post-processes.

**Why it's slow:**
- **Alloc/free per interval per iteration:** Every interval in every Baum-Welch iteration allocates fresh `alphas_1`/`betas_1` matrices. For chr1 with 10,000 intervals × 20 iterations = 200,000 heap allocations.
- **NaN/isinf checks in the hot loop:** `iForward()` and `iBackward()` contain `std::isinf`/`std::isnan` checks on every cell, which flush the FPU pipeline.
- **`long double` on x86:** `long double` is 80-bit x87 — slower than `double` on modern CPUs. Only used when `-ld` flag is set, but the default log-space code still uses it for some intermediates.

### Proposed fixes:
1. **Pre-allocate reusable buffers** — allocate one thread-local `alphas_1`/`betas_1` matrix sized to the largest interval, reuse across iterations.
2. **Guard NaN checks with compile-time debug flag** — wrap in `#ifndef NDEBUG`; they protect against a condition that shouldn't occur with correct numerics.
3. **Use `double` uniformly when `-ld` is off** — 64-bit IEEE 754 is sufficient for log-space forward-backward.
4. **Cache `logA[k_1][k_2] + eProbs[s][i][t][k_2]`** — the inner loop recomputes `logA + eProbs` for each k_1→k_2 pair; precomputing once per t saves 3× the work.

**Estimated impact:** 2–3× on forward-backward, 10–20% overall.

---

## 3. Emission probability computation — repetitive transcendental calls

**Where:** `HMM::computeEmissionProbs()` and `computeEProb()` overloads in `src/hmm_1.h` (lines ~194–440).

**What it does:** For every position t in every interval, calls:
- `gamma.getDensity(kde)` → involves `log`, `exp`, `lgamma`, and an incomplete gamma function lookup
- `bin.getDensity(count, n)` → involves `log` and `lgamma`

These are called once per position per Baum-Welch iteration. For 100,000 positions × 20 iterations = 2,000,000 calls to GSL's gamma functions.

### Proposed fixes:
1. **Use fast-math approximations** — `exp`/`log` are overkill at 80-bit precision. `__builtin_expf` or VDT (vectorized math) gives 2–4× speedup at adequate precision for likelihoods.
2. **Lazy evaluation** — `kde >= gamma1.tp` check gates the `getDensity()` call; ensure the compiler inlines this correctly at `-O3`.
3. **Move to `float` for gamma density** — the KDE values are aggregated from integer counts; `float` has sufficient dynamic range.

**Estimated impact:** 1.5–2× on emission phase, 5–10% overall.

---

## 4. GSL Nelder-Mead parameter fitting — serial + redundant

**Where:** `gamma.updateThetaAndK()` in `src/density_functions.h` (lines ~318–415).

**What it does:** Runs GSL's `gsl_multimin_fminimizer_nmsimplex2` (Nelder-Mead) to fit gamma distribution parameters to state posteriors. Each function evaluation loops over all observations. Multiple start points are tried sequentially.

**Why it's slow:**
- **Serial** — only one thread, despite being called from a parallel context.
- **Multiple start points** — tries 4 different initial (k, b0) pairs, running 4 full optimizations.
- **Full-dataset evaluations** — each Nelder-Mead step scans all observations to compute the likelihood.

### Proposed fixes:
1. **Reduce start points** — the first start point often converges to the same minimum as the others. Try 2 points instead of 4, or use the previous iteration's optimum as the single start point.
2. **Subsample for fitting** — random sample of 10–20% of observations gives essentially the same parameter estimates at a fraction of the cost.
3. **Replace Nelder-Mead with L-BFGS** — GSL's `gsl_multimin_fdfminimizer_vector_bfgs2` converges in fewer iterations when gradients are available (log-likelihood of gamma is differentiable in closed form).

**Estimated impact:** 3–5× on fitting, 5–10% overall.

---

## 5. BAM I/O — BGZF decompression

**Where:** `parse_alignments.h` via SeqAn's BAM reader.

**What it does:** Reads BAM records sequentially, builds truncCounts arrays. Currently forced to single-thread (`SEQAN_BGZF_NUM_THREADS=1` for correctness).

**Why it's slow:**
- BGZF decompression is CPU-bound
- Single-threaded access

### Proposed fix:
**Two-pass approach:** First pass reads BAM with N BGZF threads to build a raw position→count map (no HMM state, just integers). Second pass builds intervals from the map. The BGZF decompression is parallelized; the interval construction remains serial but is lightweight.

**Estimated impact:** 1.5–2× on I/O phase, 3–5% overall for large BAMs.

---

## 6. Cache efficiency — `String<String<String<double>>>` nesting

**Where:** Throughout — `eProbs`, `statePosteriors`, `setObs`, etc.

**What it does:** Uses SeqAn's `String<String<String<T>>>` for deeply nested data (strand → interval → position → state). Each `String<T>` is a heap-allocated `std::vector<T>`-like type.

**Why it's slow:**
- Four levels of pointer indirection per access
- Poor locality — interval data scattered across heap
- No vectorization possible across positions

### Proposed fix:
**Flatten to `std::vector` + offset arrays.** Replace:
```
String<String<String<double>>> statePosteriors;  // [s][k][i][t]
```
with:
```cpp
struct Interval { uint32_t offset; uint32_t length; };
std::vector<double> statePosteriors;  // flat: K*∑length
std::vector<Interval> intervals;      // [s][i]
std::vector<uint32_t> strandOffset;   // [s]
```

This provides contiguous memory access, enables SIMD, and reduces allocation count from thousands to a handful.

**Estimated impact:** 1.5–2× on HMM phases, 10–15% overall.

---

## 7. Log-sum-exp lookup table — good, but vector-hostile

**Where:** `LogSumExp_lookupTable` in `src/util.h` (lines ~42–80).

**What it does:** Precomputes `log1p(exp(x))` for 600,000 values to avoid calling `exp` + `log1p` in the inner loop. Smart approach but implemented as a single scalar lookup:

```cpp
return f1 + lookupTable[(int)(((f2-f1) - minValue) * (size/-minValue))];
```

**Why it's slow:**
- Integer-to-float conversion and array indexing in the hot path
- No SIMD — each call is scalar

### Proposed fix:
**Vectorized polynomial approximation.** Replace the lookup table with a minimax polynomial for `log1p(exp(x))` over a bounded range. This enables full SIMD vectorization and removes memory access entirely:

```cpp
// log1p(exp(x)) ≈ p(x) for x ∈ [-20, 0]
inline __m256d logSumExp_add_avx(__m256d f1, __m256d f2) {
    __m256d diff = _mm256_sub_pd(f2, f1);
    // ... polynomial eval on diff ...
    return _mm256_add_pd(f1, result);
}
```

**Estimated impact:** 2–3× on the log-sum-exp calls, 5–10% on forward-backward.

---

## Priority matrix

| # | Change | Complexity | Risk | Overall gain | Try first? |
|---|--------|-----------|------|-------------|------------|
| 1 | FFT-based KDE | Medium | Low (math-only) | **2–5×** | ✅ Yes |
| 2 | Flatten data structures | Medium | Medium (touches all code) | 10–15 % | After #1 |
| 3 | Pre-alloc FB buffers | Low | Low | 10–20 % | ✅ Yes |
| 4 | Fast-math for emission | Low | Low (guard with flag) | 5–10 % | ✅ Yes |
| 5 | Vectorized log-sum-exp | Medium | Low | 5–10 % | After #1 |
| 6 | Subsample GSL fitting | Low | Low (statistical) | 5–10 % | ✅ Yes |
| 7 | Parallel BGZF I/O | High | Medium (correctness) | 3–5 % | Last |

---

## Suggested approach

1. **Profile first** — build with `-DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_CXX_FLAGS=-pg` or use Instruments (macOS) / `perf` (Linux) on a chr21 run to validate the estimates above.
2. **Tackle the big win first** — FFT-based KDE (#1) is the single largest bottleneck and is purely a math swap; low regression risk.
3. **Then the fast wins** — pre-allocated buffers (#3), fast-math (#4), subsampling (#6) are low-complexity, low-risk, and add up to ~30% combined.
4. **Benchmark after each change** — use `ctest -L "tier1|tier2|tier3"` to verify no regression, `time` to measure wall-clock improvement.
