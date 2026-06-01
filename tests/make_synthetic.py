#!/usr/bin/env python3
"""
Generate deterministic synthetic test data for PureCLIP regression tests.

Produces:
  data/synthetic/ref.fa          – 10 000 bp single-contig reference
  data/synthetic/sample.bam      – ~880 reads with clear crosslink clusters
  data/synthetic/sample.bam.bai  – BAM index

The data is fully reproducible (seeded RNG). Re-running is a no-op if the
files already exist.

Requirements: samtools on PATH.
"""
import os, random, subprocess, sys

# ── Parameters ──────────────────────────────────────────────────────────────
SEED              = 42
CONTIG            = "chrSyn"
CONTIG_LEN        = 10_000
READ_LEN          = 25

# Forward-strand clusters: (1-based read-start position, read count)
FWD_CLUSTERS = [(500, 40), (1500, 40), (2500, 40), (3500, 40),
                (5000, 40), (6500, 40), (7500, 40), (8500, 40)]

# Reverse-strand clusters: crosslink at this position (rightmost mapped base)
# SAM POS (1-based leftmost) = xlink_pos - READ_LEN + 1
REV_CLUSTERS = [(1000, 35), (3000, 35), (5500, 35), (8000, 35)]

N_BACKGROUND      = 400   # reads spread uniformly across the contig

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "synthetic")
# ────────────────────────────────────────────────────────────────────────────


def random_seq(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGTACGTACGT") for _ in range(n))   # slight AT bias


def make_fasta(path: str, seq: str) -> None:
    with open(path, "w") as fh:
        fh.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + "\n")


def sam_read(name: str, flag: int, pos1: int, seq: str) -> str:
    """Return one SAM record line (pos1 is 1-based)."""
    qual = "I" * len(seq)
    return f"{name}\t{flag}\t{CONTIG}\t{pos1}\t255\t{len(seq)}M\t*\t0\t0\t{seq}\t{qual}"


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def make_sam(path: str, ref: str, rng: random.Random) -> None:
    records = []
    rid = 0

    # Forward-strand clusters: all reads start at the same position
    for pos1, count in FWD_CLUSTERS:
        end0 = min(pos1 - 1 + READ_LEN, CONTIG_LEN)   # 0-based exclusive
        raw  = ref[pos1 - 1 : end0].ljust(READ_LEN, "N")
        for _ in range(count):
            records.append(sam_read(f"r{rid}", 0, pos1, raw))
            rid += 1

    # Reverse-strand clusters: FLAG=16, read is reverse-complemented
    # xlink_pos (1-based) is the rightmost base of the read
    for xlink1, count in REV_CLUSTERS:
        pos1 = max(1, xlink1 - READ_LEN + 1)
        end0 = min(pos1 - 1 + READ_LEN, CONTIG_LEN)
        raw  = ref[pos1 - 1 : end0].ljust(READ_LEN, "N")
        seq  = revcomp(raw)
        for _ in range(count):
            records.append(sam_read(f"r{rid}", 16, pos1, seq))
            rid += 1

    # Background reads (forward strand)
    for _ in range(N_BACKGROUND):
        pos1 = rng.randint(1, CONTIG_LEN - READ_LEN + 1)
        end0 = pos1 - 1 + READ_LEN
        raw  = ref[pos1 - 1 : end0]
        records.append(sam_read(f"r{rid}", 0, pos1, raw))
        rid += 1

    with open(path, "w") as fh:
        fh.write(f"@HD\tVN:1.6\tSO:unsorted\n")
        fh.write(f"@SQ\tSN:{CONTIG}\tLN:{CONTIG_LEN}\n")
        for rec in records:
            fh.write(rec + "\n")


def run(cmd: str) -> None:
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"ERROR running: {cmd}\nstderr: {r.stderr}", file=sys.stderr)
        sys.exit(1)


def main() -> None:
    fa_path  = os.path.join(OUT_DIR, "ref.fa")
    sam_path = os.path.join(OUT_DIR, "sample.sam")
    bam_path = os.path.join(OUT_DIR, "sample.bam")
    bai_path = bam_path + ".bai"

    if os.path.exists(bam_path) and os.path.exists(bai_path) and os.path.exists(fa_path):
        print("Synthetic data already present — skipping generation.")
        return

    os.makedirs(OUT_DIR, exist_ok=True)

    rng = random.Random(SEED)
    ref_seq = random_seq(CONTIG_LEN, rng)

    print(f"Writing reference  → {fa_path}")
    make_fasta(fa_path, ref_seq)
    run(f"samtools faidx {fa_path}")

    print(f"Writing SAM        → {sam_path}")
    make_sam(sam_path, ref_seq, rng)

    print(f"Converting to BAM  → {bam_path}")
    run(f"samtools view -bS {sam_path} | samtools sort -o {bam_path}")
    run(f"samtools index {bam_path}")
    os.remove(sam_path)

    total_reads = sum(c for _, c in FWD_CLUSTERS) + sum(c for _, c in REV_CLUSTERS) + N_BACKGROUND
    print(f"Done — {total_reads} reads, {CONTIG_LEN} bp reference.")


if __name__ == "__main__":
    main()
