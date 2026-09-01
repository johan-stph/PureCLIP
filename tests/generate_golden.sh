#!/usr/bin/env bash
#
# Regenerate the reference ("golden") outputs that the regression tests compare
# against. Only run this when you have decided that a change to PureCLIP's
# output is correct — regenerating hides regressions.
#
# Usage:  tests/generate_golden.sh [pureclip_binary]
#
# Goldens are stored per platform ($(uname -s)) because PureCLIP's output is
# not identical across operating systems. See TESTING.md.
set -euo pipefail

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PURECLIP="${1:-${TESTS_DIR}/../build/pureclip}"
PLATFORM="$(uname -s)"

if [[ ! -x "$PURECLIP" ]]; then
    echo "ERROR: pureclip binary not found or not executable: $PURECLIP" >&2
    echo "       Build it first, or pass the path as an argument." >&2
    exit 1
fi

echo "Binary:   $PURECLIP"
echo "Platform: $PLATFORM"
echo

run_case () {
    local name="$1" data="$2" bam="$3" fa="$4"
    local out="${TESTS_DIR}/golden/${name}/${PLATFORM}"

    if [[ ! -f "${data}/${bam}" ]]; then
        echo "SKIP ${name}: ${data}/${bam} not found"
        return
    fi
    mkdir -p "$out"
    echo "=== ${name} ==="
    "$PURECLIP" \
        -i   "${data}/${bam}" \
        -bai "${data}/${bam}.bai" \
        -g   "${data}/${fa}" \
        -o   "${out}/sites.bed" \
        -or  "${out}/regions.bed" \
        -p   "${out}/params.txt" \
        -nt  4 > /dev/null
    echo "  sites:   $(wc -l < "${out}/sites.bed" | tr -d ' ') lines"
    echo "  regions: $(wc -l < "${out}/regions.bed" | tr -d ' ') lines"
}

python3 "${TESTS_DIR}/make_synthetic.py"
echo
run_case synthetic "${TESTS_DIR}/data/synthetic" sample.bam ref.fa
run_case chrM      "${TESTS_DIR}/data/chrM"      chrM.bam   chrM.fa

echo
echo "Done. Review the diff before committing:  git diff tests/golden/"
