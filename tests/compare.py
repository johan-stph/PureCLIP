#!/usr/bin/env python3
"""
PureCLIP regression test comparison utility.

Usage:
  python3 compare.py bed    <actual.bed>    <expected.bed>
  python3 compare.py params <actual.params> <expected.params> [--tol FLOAT]
  python3 compare.py count  <actual.bed>    <min_count>

Exit codes:
  0  – pass
  1  – fail
  77 – skip (golden file not found; CTest SKIP_RETURN_CODE=77)
"""
import argparse, math, os, sys


# ── BED comparison ──────────────────────────────────────────────────────────

def _sort_key(line: str):
    parts = line.split("\t")
    try:
        return (parts[0], int(parts[1]), parts[2] if len(parts) > 2 else "")
    except (IndexError, ValueError):
        return (line, 0, "")


def load_bed(path: str) -> list[str]:
    with open(path) as fh:
        return sorted(
            (l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")),
            key=_sort_key,
        )


def compare_bed(actual_path: str, expected_path: str) -> int:
    if not os.path.exists(expected_path):
        print(f"SKIP: golden BED not found: {expected_path}")
        print("      Run tests/generate_golden.sh to create golden files.")
        return 77

    actual   = load_bed(actual_path)
    expected = load_bed(expected_path)

    if actual == expected:
        print(f"OK  bed  {len(actual)} lines match.")
        return 0

    print(f"FAIL bed  actual={len(actual)} lines  expected={len(expected)} lines")
    diffs = 0
    for i, (a, e) in enumerate(zip(actual, expected)):
        if a != e:
            print(f"     line {i+1} actual:   {a}")
            print(f"     line {i+1} expected: {e}")
            diffs += 1
            if diffs >= 5:
                print(f"     … (showing first 5 differences)")
                break
    if len(actual) != len(expected):
        extra = len(actual) - len(expected)
        sign  = "extra" if extra > 0 else "missing"
        print(f"     {abs(extra)} {sign} line(s) in actual")
    return 1


# ── Params comparison ────────────────────────────────────────────────────────

def parse_params(path: str) -> dict[str, float]:
    """
    Parse a PureCLIP .params file into {key: float}.
    Lines look like:   gamma1.k<TAB>0.823
    Multi-word headers (e.g. 'Parameter learned:') are skipped.
    The transition-matrix block (tab-separated floats without a key) is
    collected as trans.0.0 … trans.3.3.
    """
    result: dict[str, float] = {}
    in_trans = False
    trans_row = 0

    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                in_trans = False
                trans_row = 0
                continue
            if line == "Transition probabilities:":
                in_trans = True
                trans_row = 0
                continue
            if in_trans:
                cols = line.split("\t")
                for col_idx, val in enumerate(cols):
                    val = val.strip()
                    if val:
                        try:
                            result[f"trans.{trans_row}.{col_idx}"] = float(val)
                        except ValueError:
                            pass
                trans_row += 1
                continue
            parts = line.split("\t")
            if len(parts) == 2:
                try:
                    result[parts[0].strip()] = float(parts[1].strip())
                except ValueError:
                    pass  # e.g. "Parameter learned:"

    return result


def compare_params(actual_path: str, expected_path: str, tol: float = 1e-4) -> int:
    if not os.path.exists(expected_path):
        print(f"SKIP: golden params not found: {expected_path}")
        print("      Run tests/generate_golden.sh to create golden files.")
        return 77

    actual   = parse_params(actual_path)
    expected = parse_params(expected_path)

    fails: list[str] = []
    all_keys = sorted(set(actual) | set(expected))

    for k in all_keys:
        if k not in actual:
            fails.append(f"  missing in actual: {k}")
        elif k not in expected:
            fails.append(f"  extra in actual:   {k}")
        else:
            a, e = actual[k], expected[k]
            denom = abs(e) if abs(e) > 1e-12 else 1.0
            rel   = abs(a - e) / denom
            if rel > tol:
                fails.append(f"  {k}: actual={a:.8g}  expected={e:.8g}  rel={rel:.2e}")

    if not fails:
        print(f"OK  params  {len(actual)} values within tol={tol}.")
        return 0

    print(f"FAIL params  {len(fails)} value(s) differ (tol={tol}):")
    for msg in fails:
        print(msg)
    return 1


# ── Count check ──────────────────────────────────────────────────────────────

def check_count(bed_path: str, min_count: int) -> int:
    with open(bed_path) as fh:
        n = sum(1 for l in fh if l.strip() and not l.startswith("#"))
    if n >= min_count:
        print(f"OK  count  {n} sites (≥ {min_count})")
        return 0
    print(f"FAIL count  {n} sites < minimum {min_count}")
    return 1


# ── CLI ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="PureCLIP regression test comparisons",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("bed", help="Compare two BED files (sorted)")
    p.add_argument("actual")
    p.add_argument("expected")

    p = sub.add_parser("params", help="Compare params files with float tolerance")
    p.add_argument("actual")
    p.add_argument("expected")
    p.add_argument("--tol", type=float, default=1e-4,
                   help="Max relative tolerance (default 1e-4)")

    p = sub.add_parser("count", help="Check BED file has at least N sites")
    p.add_argument("bed")
    p.add_argument("min_count", type=int)

    args = parser.parse_args()

    if args.cmd == "bed":
        sys.exit(compare_bed(args.actual, args.expected))
    elif args.cmd == "params":
        sys.exit(compare_params(args.actual, args.expected, args.tol))
    elif args.cmd == "count":
        sys.exit(check_count(args.bed, args.min_count))


if __name__ == "__main__":
    main()
