# Testing PureCLIP

PureCLIP fits an HMM by expectation-maximisation with a Nelder-Mead optimiser
inside it. Small numerical changes — a compiler flag, a library version, a
different floating-point width — can move the fitted parameters slightly, and
site calling amplifies those small moves into visibly different output. So the
question a test has to answer is not "did any number change" but **"does the
tool still call the same crosslink sites"**.

The reference is the current `main`. It reproduces the pre-existing chrM
reference output exactly, so it is a known-good starting point rather than an
arbitrary snapshot.

## Running the tests

```bash
make test        # synthetic fixture only, a few seconds — use during development
make test-all    # adds the real chrM data, ~30 s total — use before pushing
make test-truth  # only the ground-truth check
```

`make` builds first, so a clean checkout needs nothing else. To drive CTest
directly: `ctest --test-dir build -L tier1 --output-on-failure`.

## What the tests check

| Check | Question it answers |
|---|---|
| `smoke` | Did the run produce any output at all? |
| `sites` / `regions` | Is the output byte-identical to the reference? |
| `params` | Are the fitted model parameters within 1e-4 relative tolerance? |
| `agreement` | What fraction of called positions still agree (Jaccard ≥ 0.95)? |
| `ground_truth` | Are the implanted crosslinks recovered? |
| `determinism` | Is the output independent of `-nt` (thread count)? |

The pairing of `sites` with `agreement` is deliberate. `sites` is a tripwire: it
fails on *any* change, including one you intended. `agreement` then tells you
how bad it is — losing one site out of twelve reports `jaccard=0.9167`, whereas
a genuine breakage collapses toward zero. Read them together: `sites` says
something moved, `agreement` says whether it matters.

`ground_truth` is the only check not measured against previous output. The
synthetic fixture implants 12 crosslinks at known positions, and the test
requires all 12 to be recovered (within ±2 bp — PureCLIP reports the crosslink
one base upstream of the truncation site). It stays meaningful even if the
reference files themselves turn out to be wrong.

## Test data

- **synthetic** — 10 kb contig, 860 reads, 12 implanted crosslinks. Generated
  deterministically by `tests/make_synthetic.py` from a fixed seed, so it is
  *not* committed. Requires `samtools`.
- **chrM** — real data, committed under `tests/data/chrM/`.

## Investigating a failure

1. `make test` and read which check failed.
2. If `agreement` reports a high Jaccard, few sites moved; if it collapses,
   something structural broke.
3. If `ground_truth` fails, the change is a genuine regression in detection —
   the reference files are not involved.
4. If only `params` fails, the model fit drifted but calling survived; this is
   usually the earliest warning of a numerical change.

## Updating the reference

Only after you have confirmed the new output is correct:

```bash
make golden        # rewrites tests/golden/<platform>/
git diff tests/golden/
```

Regenerating is how a regression gets silently baptised as the new truth, so
the diff belongs in the pull request for a reviewer to look at.

## Known limitation: platform-dependent output

Goldens are stored per platform (`tests/golden/<case>/Darwin`, `.../Linux`)
because PureCLIP does not produce identical output on macOS and Linux. On chrM
the difference is large — 952 sites on macOS arm64 against 2088 on Linux — and
it is not rounding. The cause is understood:

1. `long double` differs by architecture. On macOS arm64 it *is* `double`
   (53-bit mantissa, smallest normal ~1e-308); on Linux x86-64 it is 80-bit x87
   (64-bit mantissa, ~1e-4932); on Linux aarch64 it is 128-bit quad. PureCLIP
   prints the limits at startup: on macOS arm64 you get
   `DBL_MIN_10_EXP: -307 LDBL_MIN_10_EXP: -307` — the same number twice.
2. Emission probabilities are computed in linear space and only then converted
   to logs. On a high-coverage interval the zero-truncated binomial term
   underflows to exactly 0 in `double`, while it stays representable in 80-bit.
3. When emissions hit 0, PureCLIP discards the whole interval — run with `-vv`
   to see `Warning: discarding interval ... due to emission probabilities of
   0.0`. On chrM, interval `[0, 4550)` is dropped on macOS. That window holds
   about 74% of the reads, which is why the macOS site list starts at 4893.

So macOS output is a strict *subset* of Linux output: PureCLIP under-calls
there, it never over-calls. `tests/golden/chrM/Darwin/sites.bed` is contained
exactly in the Linux list, and the same holds for regions.

PureCLIP's own advice when this happens is to rerun with `-ld` (long double).
**That flag does nothing on Apple Silicon**, because `long double` is already
`double` — the output is byte-identical and the same intervals are discarded.
The real fix is to compute the emission log-densities directly rather than
computing a linear density and taking its log; the forward-backward algorithm
is already in log space, so the emissions are the remaining gap.

Until that is fixed, compare like with like: a golden generated on macOS says
nothing about a Linux build. A missing golden makes the comparison tests
**skip** rather than fail, so a platform without recorded reference output
still gets the `ground_truth` and `determinism` coverage, which are
platform-independent.
