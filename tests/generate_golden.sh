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
CHRM_DIR="${TESTS_DIR}/data/chrM"
FULL_DIR="${TESTS_DIR}/data/full"
PLATFORM="$(uname -s)"
GOLDEN_SYN="${TESTS_DIR}/golden/synthetic/${PLATFORM}"
GOLDEN_CHR="${TESTS_DIR}/golden/chrM/${PLATFORM}"
GOLDEN_FULL="${TESTS_DIR}/golden/chr21/${PLATFORM}"

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
if [[ ! -f "${CHRM_DIR}/chrM.bam" ]]; then
    echo "SKIP: chrM data not found."
    echo "      Run:  tests/prepare_real.sh <sample_run_dir>"
else
    mkdir -p "$GOLDEN_CHR"
    "$PURECLIP" \
        -i   "${CHRM_DIR}/chrM.bam" \
        -bai "${CHRM_DIR}/chrM.bam.bai" \
        -g   "${CHRM_DIR}/chrM.fa" \
        -o   "${GOLDEN_CHR}/sites.bed" \
        -or  "${GOLDEN_CHR}/regions.bed" \
        -p   "${GOLDEN_CHR}/params.txt" \
        -nt  4

    echo "  sites:   $(wc -l < "${GOLDEN_CHR}/sites.bed")  lines"
    echo "  regions: $(wc -l < "${GOLDEN_CHR}/regions.bed") lines"
fi

# ── Tier 3: chr21 (offline, ~1 min) ─────────────────────────────────────────
echo ""
echo "=== Tier 3: chr21 golden outputs (offline) ==="
if [[ ! -f "${FULL_DIR}/chr21.bam" ]]; then
    echo "SKIP: chr21 data not found."
    echo "      Run:  tests/prepare_full.sh <sample_run_dir>"
else
    mkdir -p "$GOLDEN_FULL"
    "$PURECLIP" \
        -i   "${FULL_DIR}/chr21.bam" \
        -bai "${FULL_DIR}/chr21.bam.bai" \
        -g   "${FULL_DIR}/chr21.fa" \
        -o   "${GOLDEN_FULL}/sites.bed" \
        -or  "${GOLDEN_FULL}/regions.bed" \
        -p   "${GOLDEN_FULL}/params.txt" \
        -iv  'chr21;' -chr 'chr21;' \
        -nt  4

    echo "  sites:   $(wc -l < "${GOLDEN_FULL}/sites.bed")  lines"
    echo "  regions: $(wc -l < "${GOLDEN_FULL}/regions.bed") lines"
fi

echo ""
echo "✓ Golden outputs written."
echo "  Commit tests/golden/ to lock in expected behavior."
echo "  After optimization, 'ctest -L tier1 -L tier2' must stay green."
