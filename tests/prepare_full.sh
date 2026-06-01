#!/usr/bin/env bash
# =============================================================================
# prepare_full.sh – Extract chr21 test data for tier-3 offline validation.
#
# Run this ONCE before the tier-3 regression tests can execute.
#
# Usage:
#   ./tests/prepare_full.sh <sample_run_dir> [output_dir]
#
# Example:
#   ./tests/prepare_full.sh ../sample_run
#
# Output (in tests/data/full/ by default):
#   chr21.fa        – chr21 reference sequence
#   chr21.fa.fai    – FASTA index
#   chr21.bam       – chr21-only reads, sorted
#   chr21.bam.bai   – BAM index
# =============================================================================
set -euo pipefail

SAMPLE_RUN="${1:?Usage: $0 <sample_run_dir> [output_dir]}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${2:-${SCRIPT_DIR}/data/full}"

BAM="${SAMPLE_RUN}/aligned.prepro.bam"
REF="${SAMPLE_RUN}/ref.hg19.fa"

for f in "$BAM" "$REF"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f"; exit 1; }
done

mkdir -p "$OUT_DIR"

# ── Reference ────────────────────────────────────────────────────────────────
FA="${OUT_DIR}/chr21.fa"
if [[ ! -f "$FA" ]]; then
    echo "Extracting chr21 reference → ${FA}"
    samtools faidx "$REF" chr21 > "$FA"
    samtools faidx "$FA"
else
    echo "Reference already present: ${FA}"
fi

# ── BAM ──────────────────────────────────────────────────────────────────────
OUT_BAM="${OUT_DIR}/chr21.bam"
if [[ ! -f "$OUT_BAM" ]]; then
    if [[ ! -f "${BAM}.bai" ]]; then
        echo "Indexing source BAM (one-time, may take a minute)..."
        samtools index "$BAM"
    fi
    echo "Extracting chr21 reads → ${OUT_BAM}"
    samtools view -b "$BAM" chr21 | samtools sort -o "$OUT_BAM"
    samtools index "$OUT_BAM"
else
    echo "BAM already present: ${OUT_BAM}"
fi

echo ""
echo "✓ chr21 test data ready in ${OUT_DIR}"
echo "  Reference:  $(wc -c < "$FA" | tr -d ' ') bytes"
echo "  BAM reads:  $(samtools view -c "$OUT_BAM")"
echo ""
echo "Next step: generate tier-3 golden files with generate_golden.sh."
