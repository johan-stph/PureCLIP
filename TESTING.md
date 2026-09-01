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

Goldens are stored per platform (`tests/golden/<case>/Darwin`,
`.../Linux`) because PureCLIP does not produce identical output on macOS and
Linux. The divergence is large: on chrM the two platforms agree on the fitted
parameters to five significant figures yet call 952 versus 2088 sites. That
points at a threshold in site calling that sits close to a cliff, where a
negligible parameter difference flips many positions across the cut-off.

This is a property of PureCLIP, not of the test suite, and it is worth
investigating on its own. Until it is understood, compare like with like: a
golden generated on macOS says nothing about a Linux build. A missing golden
makes the comparison tests **skip** rather than fail, so a platform without
recorded reference output still gets the `ground_truth` and `determinism`
coverage, which are platform-independent.
