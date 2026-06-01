#!/usr/bin/env bash
# =============================================================================
# generate_golden.sh – Generate golden outputs from a reference binary.
#
# Run this ONCE (or whenever you intentionally change expected behavior) to
# create the golden files that regression tests compare against.
# Commit the resulting files in tests/golden/.
#
# Usage:
#   ./tests/generate_golden.sh <pureclip_binary> [tests_dir]
#
# Example (after building in workspace2/build):
#   ./tests/generate_golden.sh ../build/pureclip
# =============================================================================
set -euo pipefail

PURECLIP="${1:?Usage: $0 <pureclip_binary> [tests_dir]}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS_DIR="${2:-${SCRIPT_DIR}}"

SYN_DIR="${TESTS_DIR}/data/synthetic"
REAL_DIR="${TESTS_DIR}/data/real"
GOLDEN_SYN="${TESTS_DIR}/golden/synthetic"
GOLDEN_CHR="${TESTS_DIR}/golden/chrM"

[[ -x "$PURECLIP" ]] || { echo "ERROR: not executable: $PURECLIP"; exit 1; }

echo "Using binary: $PURECLIP"
echo "$($PURECLIP --version 2>&1 | head -1)"
echo ""

# ── Tier 1: Synthetic ────────────────────────────────────────────────────────
echo "=== Tier 1: synthetic golden outputs ==="
python3 "${TESTS_DIR}/make_synthetic.py"   # no-op if data exists

mkdir -p "$GOLDEN_SYN"
"$PURECLIP" \
    -i   "${SYN_DIR}/sample.bam" \
    -bai "${SYN_DIR}/sample.bam.bai" \
    -g   "${SYN_DIR}/ref.fa" \
    -o   "${GOLDEN_SYN}/sites.bed" \
    -or  "${GOLDEN_SYN}/regions.bed" \
    -p   "${GOLDEN_SYN}/params.txt" \
    -nt  4

echo "  sites:   $(wc -l < "${GOLDEN_SYN}/sites.bed")  lines"
echo "  regions: $(wc -l < "${GOLDEN_SYN}/regions.bed") lines"

# ── Tier 2: chrM ─────────────────────────────────────────────────────────────
echo ""
echo "=== Tier 2: chrM golden outputs ==="
if [[ ! -f "${REAL_DIR}/chrM.bam" ]]; then
    echo "SKIP: chrM data not found."
    echo "      Run:  tests/prepare_real.sh <sample_run_dir>"
else
    mkdir -p "$GOLDEN_CHR"
    "$PURECLIP" \
        -i   "${REAL_DIR}/chrM.bam" \
        -bai "${REAL_DIR}/chrM.bam.bai" \
        -g   "${REAL_DIR}/chrM.fa" \
        -o   "${GOLDEN_CHR}/sites.bed" \
        -or  "${GOLDEN_CHR}/regions.bed" \
        -p   "${GOLDEN_CHR}/params.txt" \
        -nt  4

    echo "  sites:   $(wc -l < "${GOLDEN_CHR}/sites.bed")  lines"
    echo "  regions: $(wc -l < "${GOLDEN_CHR}/regions.bed") lines"
fi

echo ""
echo "✓ Golden outputs written."
echo "  Commit tests/golden/ to lock in expected behavior."
echo "  After optimization, 'ctest -L tier1 -L tier2' must stay green."
