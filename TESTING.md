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

## The covariate model is compared loosely

The `-ibam` case (`tier2_cov_*`) is checked by position agreement rather than
byte equality. That model fits a GLM, and the fit amplifies small differences
in the platform's maths library: the same input on two different arm64 Macs
produces identical site positions but scores differing by 1e-4 to 1e-3
relative, and fitted parameters by up to 4e-4. The default model does not do
this — it reproduces byte-for-byte across the same two machines.

So the covariate tests assert that the same sites are called (Jaccard >= 0.98)
and that the fit is in the same place (params within 1e-2), which is what
actually has to hold. The case exists because that path was completely
untested until a crash in it shipped in v3.0.0.

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

## Cross-platform differences

Goldens are stored per platform (`tests/golden/<case>/Darwin`, `.../Linux`),
but the reason is now narrow.

PureCLIP used to build its emission probabilities in linear space and take
their logs afterwards. The zero-truncated binomial carries a `(1-p)^(n-k)`
factor that underflows to 0 once n is large, and an emission of 0 makes
PureCLIP discard the whole interval. Whether it underflowed depended on the
width of `long double` — 80-bit on x86-64 Linux, 128-bit on aarch64 Linux,
plain `double` on arm64 macOS — so macOS silently dropped high-coverage
intervals that Linux kept. On chrM it discarded `[0, 4550)`, about 74% of the
reads, and reported 952 sites where Linux reported 2088.

The emission densities are computed in log space now, so that class of
divergence is gone: macOS and Linux agree on every called position for the
chrM fixture.

What remains is small and will not go away: `lgamma`, `log` and `log1p` are not
correctly rounded, and Apple's libm differs from glibc's. Scores therefore
differ in the last digits — median relative 1.5e-5 on chrM, worst case ~2e-2.
That is enough to break a byte-exact comparison while leaving every call
identical, which is exactly why `sites` and `agreement` are paired: `sites`
will fail across platforms, `agreement` will report `jaccard=1.0000`.

So: compare like with like. A golden generated on macOS should not be used to
judge a Linux build. A missing golden makes the comparison tests **skip**
rather than fail, so a platform without recorded reference output still gets
the `ground_truth` and `determinism` coverage, which are platform-independent.
