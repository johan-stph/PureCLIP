#!/usr/bin/env python3
"""
Compare PureCLIP output against a reference, for correctness regression tests.

Usage:
  compare.py bed     <actual.bed> <expected.bed>
  compare.py agree   <actual.bed> <expected.bed> [--min-jaccard F]
  compare.py params  <actual.params> <expected.params> [--tol F]
  compare.py count   <bed> <min_count>
  compare.py recall  <sites.bed> <truth.bed> [--window N] [--min-recall F]

Exit codes: 0 pass, 1 fail, 77 skip (reference file missing).
"""
import argparse, os, sys

SKIP = 77


def _missing(path: str) -> bool:
    if not os.path.exists(path):
        print(f"SKIP  reference missing: {path}")
        return True
    return False


# ── BED helpers ─────────────────────────────────────────────────────────────
def load_bed(path: str) -> list[str]:
    """Return non-empty, non-comment lines with trailing whitespace stripped."""
    with open(path) as fh:
        return [ln.rstrip() for ln in fh
                if ln.strip() and not ln.startswith(("#", "track", "browser"))]


def sites(path: str) -> set[tuple[str, str, str]]:
    """Reduce a BED to a set of (contig, start, strand) site keys.

    Score and name are deliberately ignored: this set answers "were the same
    positions called", which is the question that survives numerical drift.
    Strand is column 6 when present.
    """
    out = set()
    for ln in load_bed(path):
        f = ln.split("\t")
        if len(f) < 2:
            continue
        out.add((f[0], f[1], f[5] if len(f) > 5 else "."))
    return out


# ── Modes ───────────────────────────────────────────────────────────────────
def compare_bed(actual_path: str, expected_path: str) -> int:
    """Exact comparison. The strict gate: any change at all is reported."""
    if _missing(expected_path):
        return SKIP
    actual, expected = load_bed(actual_path), load_bed(expected_path)
    if actual == expected:
        print(f"OK    bed     {len(actual)} line(s) identical.")
        return 0

    a_set, e_set = set(actual), set(expected)
    print(f"FAIL  bed     {os.path.basename(actual_path)} differs from reference:")
    print(f"        actual {len(actual)} line(s), expected {len(expected)}")
    only_a, only_e = a_set - e_set, e_set - a_set
    for label, rows in (("only in actual", only_a), ("only in expected", only_e)):
        if rows:
            print(f"        {len(rows)} {label}:")
            for r in sorted(rows)[:5]:
                print(f"          {r}")
            if len(rows) > 5:
                print(f"          ... and {len(rows) - 5} more")
    if not only_a and not only_e:
        print("        same lines, different order")
    return 1


def compare_agree(actual_path: str, expected_path: str, min_jaccard: float) -> int:
    """Set agreement on called positions, tolerant of score/ordering drift."""
    if _missing(expected_path):
        return SKIP
    a, e = sites(actual_path), sites(expected_path)
    if not a and not e:
        print("OK    agree   both empty.")
        return 0

    inter, union = len(a & e), len(a | e)
    jaccard = inter / union if union else 1.0
    recall = inter / len(e) if e else 1.0
    precision = inter / len(a) if a else 1.0
    print(f"      agree   jaccard={jaccard:.4f} recall={recall:.4f} "
          f"precision={precision:.4f} (actual={len(a)} reference={len(e)})")
    if jaccard < min_jaccard:
        print(f"FAIL  agree   jaccard {jaccard:.4f} < required {min_jaccard}")
        return 1
    print(f"OK    agree   jaccard >= {min_jaccard}")
    return 0


def parse_params(path: str) -> dict[str, float]:
    """Read key/value params; unnamed rows (the transition matrix) get indexed keys."""
    values, row = {}, 0
    with open(path) as fh:
        for ln in fh:
            f = ln.strip().split("\t")
            f = [x for x in f if x]
            if not f:
                continue
            if len(f) == 2:
                try:
                    values[f[0]] = float(f[1])
                    continue
                except ValueError:
                    pass
            try:
                floats = [float(x) for x in f]
            except ValueError:
                continue
            for col, v in enumerate(floats):
                values[f"matrix[{row}][{col}]"] = v
            row += 1
    return values


def compare_params(actual_path: str, expected_path: str, tol: float) -> int:
    """Relative-tolerance comparison of the fitted model parameters."""
    if _missing(expected_path):
        return SKIP
    actual, expected = parse_params(actual_path), parse_params(expected_path)

    missing = set(expected) - set(actual)
    added = set(actual) - set(expected)
    if missing or added:
        print("FAIL  params  key mismatch")
        if missing:
            print(f"        missing: {sorted(missing)[:8]}")
        if added:
            print(f"        unexpected: {sorted(added)[:8]}")
        return 1

    fails, worst = [], 0.0
    for k, e in expected.items():
        a = actual[k]
        denom = abs(e) if abs(e) > 1e-12 else 1.0
        rel = abs(a - e) / denom
        worst = max(worst, rel)
        if rel > tol:
            fails.append(f"  {k}: actual={a:.8g} expected={e:.8g} rel={rel:.2e}")

    if fails:
        print(f"FAIL  params  {len(fails)} of {len(expected)} value(s) exceed tol={tol}:")
        for line in fails[:10]:
            print(line)
        if len(fails) > 10:
            print(f"  ... and {len(fails) - 10} more")
        return 1
    print(f"OK    params  {len(actual)} values within tol={tol} (worst rel={worst:.2e}).")
    return 0


def check_count(bed_path: str, min_count: int) -> int:
    """Smoke check: the run produced a plausible number of calls."""
    n = len(load_bed(bed_path))
    if n < min_count:
        print(f"FAIL  count   {os.path.basename(bed_path)} has {n} line(s), expected >= {min_count}")
        return 1
    print(f"OK    count   {n} line(s) >= {min_count}.")
    return 0


def check_recall(sites_path: str, truth_path: str, window: int, min_recall: float) -> int:
    """Score calls against implanted crosslinks (synthetic data only).

    Unlike every other mode this does not compare against previous output, so
    it stays meaningful even when the reference itself is wrong.
    """
    if _missing(truth_path):
        return SKIP
    called = [(c, int(s), st) for c, s, st in sites(sites_path) if s.isdigit()]
    truth = [(c, int(s), st) for c, s, st in sites(truth_path) if s.isdigit()]
    if not truth:
        print(f"FAIL  recall  no truth records in {truth_path}")
        return 1

    hits = sum(
        any(c == tc and st == tst and abs(s - ts) <= window for c, s, st in called)
        for tc, ts, tst in truth
    )
    recall = hits / len(truth)
    print(f"      recall  {hits}/{len(truth)} implanted crosslinks recovered "
          f"(window=+/-{window}bp, {len(called)} sites called)")
    if recall < min_recall:
        print(f"FAIL  recall  {recall:.4f} < required {min_recall}")
        return 1
    print(f"OK    recall  {recall:.4f} >= {min_recall}")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("bed");    p.add_argument("actual"); p.add_argument("expected")
    p = sub.add_parser("agree");  p.add_argument("actual"); p.add_argument("expected")
    p.add_argument("--min-jaccard", type=float, default=0.95)
    p = sub.add_parser("params"); p.add_argument("actual"); p.add_argument("expected")
    p.add_argument("--tol", type=float, default=1e-4)
    p = sub.add_parser("count");  p.add_argument("bed"); p.add_argument("min_count", type=int)
    p = sub.add_parser("recall"); p.add_argument("sites"); p.add_argument("truth")
    p.add_argument("--window", type=int, default=2)
    p.add_argument("--min-recall", type=float, default=1.0)

    a = parser.parse_args()
    if a.cmd == "bed":
        sys.exit(compare_bed(a.actual, a.expected))
    if a.cmd == "agree":
        sys.exit(compare_agree(a.actual, a.expected, a.min_jaccard))
    if a.cmd == "params":
        sys.exit(compare_params(a.actual, a.expected, a.tol))
    if a.cmd == "count":
        sys.exit(check_count(a.bed, a.min_count))
    if a.cmd == "recall":
        sys.exit(check_recall(a.sites, a.truth, a.window, a.min_recall))


if __name__ == "__main__":
    main()
