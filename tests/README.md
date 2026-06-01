# PureCLIP Regression Test Suite

30 tests across 3 tiers. Tier-1 + tier-2 run in CI on every push/PR.
Tier-3 is offline-only (needs external sample data).

## Overview

| Tier | Data | Tests | Wall time | RAM | CI |
|------|------|-------|-----------|-----|----|
| 1 | Synthetic 10 kb, 860 reads | 15 | ~4 s | < 256 MB | ✅ macOS + Linux |
| 2 | Real chrM 16 kb, 97k reads (committed, 3 MB) | 7 | ~12 s | < 512 MB | ✅ macOS + Linux |
| 3 | Real chr21 48 Mb, 77k reads (offline, 49 MB ref) | 6 | ~50 s | < 1 GB | ❌ manual only |

---

## Quick start — CI-level tests (no setup)

```bash
cd PureCLIP-workspace2
mkdir -p build && cd build
cmake ../src -DCMAKE_BUILD_TYPE=Release
make -j$(sysctl -n hw.logicalcpu)
ctest -L "tier1|tier2" --output-on-failure -j4
```

Expect: `100% tests passed, 0 tests failed out of 22` in ~12 seconds.

---

## Enabling tier-3 (one-time setup)

```bash
# Extract chr21 from your sample data (49 MB reference, 77k reads)
cd PureCLIP-workspace2
./tests/prepare_full.sh ../sample_run

# (Re)generate all golden files from the reference binary
./tests/generate_golden.sh build/pureclip

# Re-configure to pick up the new data
cd build && cmake ../src
# → Tier 3 tests ENABLED (chr21 data found)
```

---

## Developer workflow — real examples

### After every build: smoke test
```bash
ctest -L tier1 -j4 --output-on-failure
# 100% tests passed, 0 tests failed out of 15   ← ~4 s
```

### Before pushing a PR: full CI suite
```bash
ctest -L "tier1|tier2" -j4 --output-on-failure
# 100% tests passed, 0 tests failed out of 22   ← ~12 s
```

### Before merging: complete pre-merge validation
```bash
ctest -L "tier1|tier2|tier3" -j4 --output-on-failure
# 100% tests passed, 0 tests failed out of 28   ← ~50 s
```

### Fast regressions only (skip smoke/mode tests)
```bash
ctest -L regression -j4 --output-on-failure
# 100% tests passed, 0 tests failed out of 9    ← ~50 s
# Instantly catches any output-drifting change.
```

### Run a single failing test with full output
```bash
ctest -R t1_2_regression_sites --output-on-failure
# FAIL  bed  actual=12 lines  expected=12 lines
#      line 3 actual:   chrSyn  1498  1499  3  90.9891  ...
#      line 3 expected: chrSyn  1498  1499  3  90.8621  ...
```

### Thread determinism check
```bash
ctest -L determinism --output-on-failure
# 100% tests passed, 0 tests failed out of 1    ← ~1 s
# Confirms -nt 4 produces bit-identical BED as -nt 1.
```

---

## What happens when a test fails

All regression tests do **exact BED comparison** (lines sorted by chr:start).
Even a 0.1% shift in crosslink scores gets caught:

```
FAIL bed  actual=12 lines  expected=12 lines
     line 3 actual:   chrSyn  1498  1499  3  90.9891  +   [score_CL=250.972;score_E=90.9891;score_B=341.961;score_UC=90.9891]
     line 3 expected: chrSyn  1498  1499  3  90.8621  +   [score_CL=251.062;score_E=90.8621;score_B=341.924;score_UC=90.8621]
```

Parameter regression uses float tolerance (`tol=1e-4` by default):

```
FAIL params  3 value(s) differ (tol=0.0001):
  gamma2.b0: actual=-1.86378  expected=-1.87234  rel=4.57e-03
  gamma2.mean: actual=0.155085  expected=0.153764  rel=8.59e-03
  gamma2.k: actual=2.04298  expected=1.99532  rel=2.39e-02
```

---

## Golden files

Golden files live in `tests/golden/` and are committed to the repo.
They are the ground-truth outputs produced by the **reference binary**
(before any optimisation).

**Regenerate after intentional output changes:**
```bash
./tests/generate_golden.sh build/pureclip
git add tests/golden/
git commit -m "chore: update golden outputs"
```

**Missing golden files cause tests to SKIP (exit 77), not FAIL.**

---

## Test list

### Tier 1 — Synthetic (15 tests, ~4 s)

| # | Name | Assertion |
|----|------|-----------|
| — | `t1_0_make_data` | Generate synthetic BAM/FASTA (setup) |
| — | `t1_1_smoke` + `_check` | Exits 0, sites.bed ≥ 1 line |
| ⚡ | `t1_2_regression_sites` | sites.bed **exact match** vs golden |
| ⚡ | `t1_3_regression_regions` | regions.bed **exact match** vs golden |
| ⚡ | `t1_4_regression_params` | params.txt float match within `1e-4` |
| — | `t1_5_st1` / `st2` / `st3` | Scoring type 1/2/3 exits 0 |
| — | `t1_6_bc1` | Binding characteristics mode 1 exits 0 |
| — | `t1_7_outputall` + `_check` | `-oa` flag, sites ≥ 1 |
| 🔀 | `t1_8_nt4` + `_thread_determinism` | `-nt 4` BED ≡ `-nt 1` golden |

⚡ = regression  🔀 = determinism

### Tier 2 — Real chrM (7 tests, ~12 s, committed data)

| # | Name | Assertion |
|----|------|-----------|
| — | `t2_1_smoke_chrM` + `_check` | Real data, exits 0, sites ≥ 1 |
| ⚡ | `t2_2_regression_chrM_sites` | sites.bed **exact match** vs golden |
| ⚡ | `t2_3_regression_chrM_regions` | regions.bed **exact match** vs golden |
| ⚡ | `t2_4_regression_chrM_params` | params.txt float match within `1e-4` |
| — | `t2_5_chrM_bc1` | bc=1 mode on real data exits 0 |

### Tier 3 — Real chr21 (6 tests, ~50 s, needs `prepare_full.sh`)

| # | Name | Assertion |
|----|------|-----------|
| — | `t3_1_smoke_chr21` + `_check` | 48 Mb chromosome, exits 0, sites ≥ 1 |
| ⚡ | `t3_2_regression_chr21_sites` | sites.bed **exact match** vs golden |
| ⚡ | `t3_3_regression_chr21_regions` | regions.bed **exact match** vs golden |
| ⚡ | `t3_4_regression_chr21_params` | params.txt float match within `1e-4` |

---

## All CTest invocations

```bash
ctest -N                                  # dry-run: list all 30 tests
ctest -L tier1                            # fast sanity check (~4 s)
ctest -L "tier1|tier2"                    # CI suite (~12 s)
ctest -L "tier1|tier2|tier3"              # full pre-merge validation (~50 s)
ctest -L regression                       # golden-comparison only
ctest -L determinism                      # thread consistency
ctest -R t1_2_regression_sites            # single failing test
ctest --output-on-failure -j4             # parallel, show diffs on fail
ctest --repeat until-fail:3               # flake catcher
```

---

## compare.py (standalone)

```bash
python3 tests/compare.py bed    actual.bed    golden/sites.bed
python3 tests/compare.py params actual.params golden/params.txt --tol 1e-4
python3 tests/compare.py count  actual.bed    10
```

Exit codes: `0` = pass, `1` = fail, `77` = skip (golden missing → CTest `SKIP_RETURN_CODE`).
