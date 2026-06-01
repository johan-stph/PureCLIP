#!/usr/bin/env bash
# =============================================================================
# prepare_real.sh – Extract chrM test data from a sample_run directory.
#
# Run this ONCE manually before the tier-2 regression tests can execute.
#
# Usage:
#   ./tests/prepare_real.sh <sample_run_dir> [output_dir]
#
# Example:
#   ./tests/prepare_real.sh ../sample_run
#
# Output (in tests/data/real/ by default):
#   chrM.fa        – chrM reference sequence
#   chrM.fa.fai    – FASTA index
#   chrM.bam       – chrM-only reads, sorted
#   chrM.bam.bai   – BAM index
# =============================================================================
set -euo pipefail

SAMPLE_RUN="${1:?Usage: $0 <sample_run_dir> [output_dir]}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${2:-${SCRIPT_DIR}/data/chrM}"

BAM="${SAMPLE_RUN}/aligned.prepro.bam"
REF="${SAMPLE_RUN}/ref.hg19.fa"

# ── Validate inputs ──────────────────────────────────────────────────────────
for f in "$BAM" "$REF"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f"; exit 1; }
done

mkdir -p "$OUT_DIR"

# ── Reference ────────────────────────────────────────────────────────────────
FA="${OUT_DIR}/chrM.fa"
if [[ ! -f "$FA" ]]; then
    echo "Extracting chrM reference → ${FA}"
    samtools faidx "$REF" chrM > "$FA"
    samtools faidx "$FA"
else
    echo "Reference already present: ${FA}"
fi

# ── BAM ──────────────────────────────────────────────────────────────────────
OUT_BAM="${OUT_DIR}/chrM.bam"
if [[ ! -f "$OUT_BAM" ]]; then
    # Index the source BAM if needed (slow-path idxstats already walked it)
    if [[ ! -f "${BAM}.bai" ]]; then
        echo "Indexing source BAM (one-time, may take a minute)..."
        samtools index "$BAM"
    fi
    echo "Extracting chrM reads → ${OUT_BAM}"
    samtools view -b "$BAM" chrM | samtools sort -o "$OUT_BAM"
    samtools index "$OUT_BAM"
else
    echo "BAM already present: ${OUT_BAM}"
fi

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "✓ chrM test data ready in ${OUT_DIR}"
echo "  Reference:  $(wc -c < "$FA" | tr -d ' ') bytes"
echo "  BAM reads:  $(samtools view -c "$OUT_BAM")"
echo ""
echo "Next step: run generate_golden.sh to create the golden outputs."
