# PureCLIP Regression Test Suite

## Overview

| Tier | Data | Tests | Typical time | RAM |
|------|------|-------|-------------|-----|
| 1 | Synthetic (10 kb, ~880 reads) | 13 | ~1–5 s each | < 256 MB |
| 2 | Real chrM (16 kb, ~97 k reads) | 6 | ~30–60 s each | < 512 MB |

Tier-2 tests are automatically disabled at cmake-configure time if
`tests/data/real/chrM.bam` does not exist.

## Quick Start

```bash
# 1. Build
cd PureCLIP-workspace2
mkdir build && cd build
cmake ../src
make -j$(sysctl -n hw.logicalcpu)

# 2. Run tier-1 tests (fast, no data prep needed)
ctest -L tier1 --output-on-failure

# 3. (Optional) Enable tier-2 real-data tests
cd ..
./tests/prepare_real.sh ../sample_run    # extract chrM once

# 4. Re-configure so CMake sees the new data, then test
cd build && cmake ../src
ctest -L "tier1|tier2" --output-on-failure
```

## Golden Files

Golden files live in `tests/golden/` and are committed to the repo.
They are the ground-truth outputs produced by the **reference binary**
(before any optimisation).

Generate or re-generate them with:
```bash
./tests/generate_golden.sh build/pureclip
git add tests/golden/
git commit -m "chore: update golden outputs"
```

**Missing golden files cause tests to SKIP (not FAIL).** This lets you run
smoke and mode tests on a fresh checkout without needing to regenerate.

## Test List

### Tier 1 – Synthetic

| Test | What it checks |
|------|----------------|
| `t1_0_make_data` | Generates `data/synthetic/` (setup) |
| `t1_1_smoke` | Basic run exits 0 and produces non-empty output |
| `t1_1_smoke_check` | sites.bed has ≥ 1 line |
| `t1_2_reg_run` | Regression run (reused by t1_2/3/4) |
| `t1_2_regression_sites` | sites.bed **exact** match vs golden (sorted) |
| `t1_3_regression_regions` | regions.bed **exact** match vs golden |
| `t1_4_regression_params` | params.txt float values within `tol=1e-4` |
| `t1_5_st1/2/3` | Scoring type variants exit 0 |
| `t1_6_bc1` | Binding-characteristics mode 1 exits 0 |
| `t1_7_outputall` | `-oa` flag exits 0 and produces ≥ 1 site |
| `t1_8_nt4` + `t1_8_thread_determinism` | 4-thread run matches golden sites.bed |

### Tier 2 – Real chrM

| Test | What it checks |
|------|----------------|
| `t2_1_smoke_chrM` | Basic run on real data exits 0 |
| `t2_2_reg_run` | Regression run (reused by t2_2/3/4) |
| `t2_2_regression_chrM_sites` | sites.bed **exact** match vs golden |
| `t2_3_regression_chrM_regions` | regions.bed **exact** match vs golden |
| `t2_4_regression_chrM_params` | params float values within `tol=1e-4` |
| `t2_5_chrM_bc1` | bc=1 mode on real data exits 0 |

## Useful CTest Invocations

```bash
ctest -L tier1                            # fast checks only
ctest -L "tier1|tier2"                    # full suite
ctest -L regression                       # only golden comparisons
ctest -L determinism                      # thread consistency
ctest -N                                  # dry-run: list all tests
ctest --output-on-failure -j4             # parallel, show failures
ctest -R t1_2                             # run tests matching regex
```

## compare.py

The comparison utility is standalone and can be called directly:

```bash
python3 tests/compare.py bed    actual.bed    golden/synthetic/sites.bed
python3 tests/compare.py params actual.params golden/synthetic/params.txt --tol 1e-4
python3 tests/compare.py count  actual.bed    10
```

Exit codes: `0` = pass, `1` = fail, `77` = skip (golden missing).
