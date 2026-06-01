# PureCLIP Performance Improvement Opportunities

Analysis of `src/` — ~7,000 lines of C++ header-only HMM pipeline.
Benchmark target: chr21 (49 MB reference, 77k reads) on macOS ARM64 (Apple M-series).
Profiling tool: `sample` (macOS time profiler).

---

## ✅ Resolved — chr21 tier3 wall time: 86.2s → 48.0s (−44.3%)

All optimizations verified with `ctest -L "tier1|tier2|tier3"` (21 tests, all pass).
Output is ULPs-identical to baseline (86 sites, 42 regions; parameter drift at 4th–5th decimal).

### 1. KDE convolution — SIMD + FFT (was #1)

**Commit:** `4f7951a`, `1eba0af`  
**Impact:** 86.2s → 73.3s (−15.1%)

- Split the O(T×W) windowed sum into two branchless halves with `#pragma omp simd reduction`  
- Added vDSP FFT path (Accelerate framework) for intervals > 500 bp on macOS  
- FFT was marginal on chr21 (intervals 500–4000 bp, 7300 FFT invocations) — overhead dominates. Would shine on chr1 (249 Mb) where intervals are 10k–50k bp.

### 2. Forward-backward — pre-allocated buffers + NDEBUG guards (was #3)

**Commit:** `dffbabf`  
**Impact:** 52.5s → 49.6s (−5.7%)

- Pre-allocate thread-local `alphas_1`/`betas_1`/`xis`/`p_i` once per parallel region (sized to max interval), eliminating ~200k heap allocs per run  
- Wrapped `std::isinf`/`std::isnan` checks in `#ifndef NDEBUG` — these flush the FPU pipeline on every cell

### 3. GSL Nelder-Mead fitting — fewer start points + relaxed tolerance (was #4, #6)

**Commit:** `27eeccb`, `805f9af`  
**Impact:** 73.3s → 52.5s (−28.4%)

- Reduced start points from 4–6 → 2 (saves 50% of simplex runs; converges to same minimum)  
- Reduced `maxIter_simplex` from 2000 → 500 (2-param Nelder-Mead converges much sooner)  
- Relaxed simplex size tolerance from 1e-6 → 1e-4 (coarser convergence during intermediate Baum-Welch steps)  
- Subsampling infrastructure is in place (`bool subsample` parameter) but disabled — 33% subsample lost 3 crosslink sites

### 4. Emission probs — gamma density cache (was #4)

**Commit:** `480eab8`  
**Impact:** 49.6s → 48.0s (−3.1%)

- Cached `theta`, `pow(theta,k)`, `tgamma(k)`, `gamma_p(k, tp/theta)` in `GAMMA` — recomputed lazily when `b0`/`k` change (once per Baum-Welch parameter update, not per position)  
- Per-position cost drops from 2×pow + tgamma + gamma_p + exp → 1×exp + 1×log

---

## 🔬 Measured Profile — chr21 After All Resolved Optimizations

15,331 samples on main thread, 20s window during `learnModel` (macOS `sample`):

```
learnModel() ──────────────────────────── 100%
├── learnHMM #1 ────────────────────────  68%
│   ├── baumWelch ──────────────────────  44%
│   │   ├── GSL simplex (gamma fit) ───  25%
│   │   └── forward-backward ──────────  16%
│   └── computeEmissionProbs ───────────   8%
├── learnHMM #2 ────────────────────────  32%
│   ├── computeEmissionProbs ───────────  10%
│   └── baumWelch ──────────────────────  rest
└── prior_mle (GSL simplex init) ────────  11%
```

**Hot leaves (all threads):**

| Function | Samples | Category |
|----------|---------|----------|
| `my_GSL_X_GAMMA_forK` | 21,823 | Gamma fitting (already optimized 50%) |
| `computeStatePosteriorsFBupdateTrans` | 6,730 | Forward-backward (pre-alloc done, NDEBUG done) |
| `__psynch_cvwait` (thread idle) | 4,822 | **Sync overhead** |
| `exp` | 4,102 | **Emission + FB (myExp)** |
| `log` | 2,478 | **Emission + FB (myLog)** |
| `pow` | 2,012 | **Emission (ZTBIN binomial)** |
| `tgamma` | 1,685 | Gamma density (already cached) |
| `_xzm_free` (malloc) | 1,560 | **Alloc churn** |
| `Fct_GSL_X_GAMMA_fixK` | 1,449 | GSL fixK simplex |
| `iBackward` + `iForward` | 3,111 | FB log-sum-exp (already NDEBUG'd) |
| `get_logSumExp_states` | 260 | Log-sum-exp lookup table |
| `computeEProb` | 576 | Emission probability dispatch |
| `GAMMA::getDensity` | 488 | Gamma PDF (already cached) |

**→ Remaining top targets: ZTBIN binomial (pow/exp), thread sync (cvwait), alloc churn**

---

## 📋 TODO — Remaining Opportunities (ordered by estimated impact)

### A. ZTBIN binomial density — rebuild per-call (emission hot path)

**Where:** `ZTBIN::getDensity()` in `src/density_functions_crosslink.h:95`

**Profile evidence:** `pow` at 2,012 samples + `boost::math::pdf` underlying `exp`. Called 2× per position.

**What it does:** Constructs a new `boost::math::binomial_distribution<long double>` object on every call, then calls `pdf()`, then computes `pow(1-p, n)` for zero-truncation correction.

**Proposed fixes:**
1. **Cache binomial distribution per `n` value** — `nEstimates` has limited cardinality (typically 1–50). Store a map/pool of pre-built distributions keyed by `n`, rebuild only when `p` changes.
2. **Replace `pow(1-p, n)` with `exp(n * log1p(-p))`** — faster and more accurate for small `p`.
3. **Use `double` instead of `long double`** — 64-bit is sufficient for binomial likelihood.

**Estimated impact:** 2–4% overall.

---

### B. Log-sum-exp lookup — vectorize or replace with polynomial

**Where:** `LogSumExp_lookupTable` / `get_logSumExp_states()` in `src/util.h:42`

**Profile evidence:** 260 samples in `get_logSumExp_states` itself, but called from inner loops of iForward/iBackward (3,111 combined samples).

**What it does:** Precomputed `log1p(exp(x))` table with 600,000 entries. Scalar integer-to-float conversion + array indexing per call.

**Proposed fix:** Replace with a **minimax polynomial approximation** for `log1p(exp(x))` over x ∈ [−20, 0]. Enables full SIMD vectorization, eliminates memory access.

```cpp
// Fast approximate log1p(exp(x)) via polynomial (max error ~1e-7)
inline double fast_logsumexp_add(double x) {
    if (x > 0) return x + log1p(exp(-x));  // large x: use identity
    // Polynomial approximation for x in [-20, 0]
    const double c[] = {...};  // minimax coefficients
    double x2 = x * x;
    return ((c[0] * x2 + c[1]) * x2 + c[2]) * x2 + c[3];
}
```

**Estimated impact:** 3–5% overall (accelerates both FB and emission hot paths).

---

### C. Thread synchronization overhead

**Profile evidence:** `__psynch_cvwait` at 4,822 samples (14% of main thread time).

**What it does:** OpenMP threads spin-waiting at implicit barriers — `parallel for` with `schedule(dynamic, 1)` creates a barrier at loop end. With varying interval sizes, threads finish at different times.

**Proposed fixes:**
1. **Use `schedule(guided)` instead of `schedule(dynamic, 1)`** — larger initial chunks reduce barrier frequency.
2. **Merge adjacent parallel regions** — combine `computeEmissionProbs` + `computeStatePosteriorsFBupdateTrans` into one parallel region to halve barrier count.
3. **Use `nowait` on non-dependent loops** — some loops can proceed without waiting for all threads.

**Estimated impact:** 5–10% overall.

---

### D. ZTBIN binomial — use `double` instead of `long double`

**Where:** `ZTBIN::getDensity()` and `ZTBIN::updateP()` in `src/density_functions_crosslink.h`

**Profile evidence:** `pow` at 2,012 samples (mostly from `pow(1-p, n)` in binomial). ARM64 has native 64-bit FP; 80-bit `long double` is emulated in software.

**Proposed fix:** Switch binomial computations to `double` — the input values (counts, n estimates) are integers in range 0–255, so 64-bit precision is more than sufficient.

**Estimated impact:** 2–5% overall.

---

### E. Memory allocation churn

**Profile evidence:** `_xzm_free` at 1,560 samples + `_xzm_xzone_malloc` at 745 + `free` at 533 = 2,838 samples (~8%).

**What it does:** Per-interval allocations for `eProbs`, `statePosteriors`, `initProbs` arrays. These are `String<String<String<double>>>` nested types — each inner `String<T>` is a separate heap allocation.

**Proposed fix:** Flatten `eProbs` to contiguous `std::vector<double>` with offset arrays (as described in original #6). Also pre-allocate emission probability buffers per thread like we did for alphas/betas.

**Estimated impact:** 5–10% overall.

---

### F. BAM I/O — parallel BGZF decompression

**Where:** `parse_alignments.h`, forced `SEQAN_BGZF_NUM_THREADS=1`

**Profile evidence:** Not profiled on chr21 (I/O is <5% of 48s runtime). Would matter on full-genome runs.

**Proposed fix:** Two-pass approach — first pass reads BAM with N threads to build position→count map, second pass builds intervals serially.

**Estimated impact:** 1–3% on chr21, 5–10% on full-genome.

---

### G. Full-genome scaling — vDSP FFT for large intervals

**Profile evidence:** FFT path already implemented but chr21 intervals are too small (500–4000 bp) for the O(T log T) advantage to overcome overhead. On chr1 (249 Mb), intervals would be 10k–50k bp where FFT dominates.

**Proposed fix:** Already implemented (`1eba0af`). Lower the threshold from 500→200 on production runs with large chromosomes. Pre-allocate FFT buffers (currently alloc per call).

**Estimated impact:** Negligible on chr21; 10–30% on full-genome.

---

### H. Replace `long double` with `double` throughout (when `-ld` is off)

**Where:** Throughout — FB, emission, density functions.

**Profile evidence:** ARM64 has no native 80-bit FPU; `long double` operations are software-emulated. The current code uses `long double` extensively even when the `-ld` flag (high precision mode) is off.

**Proposed fix:** Use `double` for all computation when not in high-precision mode. Keep `long double` only behind `#ifdef USE_LONG_DOUBLE`.

**Estimated impact:** 10–20% overall (ARM64-specific; negligible on x86-64 where 80-bit x87 is hardware).

---

## Updated Priority Matrix

| # | Change | Complexity | Risk | Est. gain | Try? |
|---|--------|-----------|------|-----------|------|
| A | Cache ZTBIN binomial dist | Low | Low | 2–4% | ✅ |
| B | Vectorized log-sum-exp poly | Medium | Low | 3–5% | ✅ |
| C | Merge parallel regions + guided schedule | Low | Low | 5–10% | ✅ |
| D | ZTBIN double instead of long double | Low | Low | 2–5% | ✅ |
| E | Flatten allocations (eProbs/statePost) | Medium | Medium | 5–10% | After A–D |
| F | Parallel BGZF I/O | High | Medium | 1–3% (chr21) | Last |
| G | vDSP FFT for large intervals | Done | — | 0% (chr21) | Ready for production |
| H | Double instead of long double | High | Medium | 10–20% | Investigate |

---

## Suggested Next Steps

1. **Tackle A+B+C+D** — all low-risk, low-complexity, combined estimate 12–24% additional reduction  
2. **Profile after A–D** — see if alloc churn (E) is now the top bottleneck  
3. **H (double vs long double)** — if profiling confirms ARM64 software emulation overhead, this is the single biggest remaining win  
4. **Full-genome benchmark** — validate G (FFT) and F (BAM I/O) on chr1 or whole-genome data
