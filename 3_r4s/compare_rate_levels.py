#!/usr/bin/env python3
"""Cross-level comparison of Rate4Site runs (Chordopoxvirinae vs Orthopoxvirus).

Rate4Site normalises its scores within each run (mean 0, sd 1 over the reference
positions), so absolute scores are not comparable between runs. Only the relative
ordering of positions is. This script therefore reports:

  * Spearman rank correlation between the two runs, which is invariant to the
    normalisation and is the primary cross-level statistic;
  * a paired Wilcoxon signed-rank test on the within-run ranks, which asks whether
    positions systematically change their relative position between levels;
  * the same test on the normalised scores, for completeness.

Positions flagged low-confidence (wide 25-75% posterior interval) can be excluded
with --drop-flagged; the flagged set is read from the .res files themselves.

Usage:
    python3 compare_rate_levels.py \
        --a r4s_run_subfamily/r60_jtt.res  --label-a "60 seq (Chordopoxvirinae)" \
        --b r4s_run_orthopox/r14_jtt.res   --label-b "14 seq (Orthopoxvirus)" \
        [--by-domain] [--drop-flagged] [--csv rate_comparison.csv]
"""

import argparse
import math
import re
import sys

DOMAINS = [
    ("N-terminal region",      1,  20),
    ("GAG-binding domain",     21, 32),
    ("Spacer",                 33, 42),
    ("Coiled-coil domain",     43, 84),
    ("Disulfide-linked motif", 71, 72),
    ("Leucine zipper domain",  85, 110),
]

ROW = re.compile(
    r"^\s*(\d+)\s+(\S)\s+(-?[\d.]+(?:[eE][-+]?\d+)?)\s*"
    r"(?:\[\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)\s*\])?"
)


def read_res(path):
    """Parse a Rate4Site .res file -> {pos: dict(aa, score, lo, hi)}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            m = ROW.match(line)
            if not m:
                continue
            pos, aa, score, lo, hi = m.groups()
            out[int(pos)] = {
                "aa": aa,
                "score": float(score),
                "lo": float(lo) if lo is not None else None,
                "hi": float(hi) if hi is not None else None,
            }
    if not out:
        sys.exit(f"no data rows parsed from {path}")
    return out


def flagged_positions(res, fraction=1 / 3):
    """ConSurf-style low-confidence flag: QQ interval wider than a third of range."""
    scores = [r["score"] for r in res.values()]
    span = max(scores) - min(scores)
    flagged = set()
    for pos, r in res.items():
        if r["lo"] is None or r["hi"] is None:
            continue
        if (r["hi"] - r["lo"]) > fraction * span:
            flagged.add(pos)
    return flagged


def rankdata(values):
    """Average ranks, ties shared."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    if den == 0:
        return float("nan"), float("nan")
    rho = num / den
    if n < 4 or abs(rho) >= 1:
        return rho, float("nan")
    t = rho * math.sqrt((n - 2) / (1 - rho ** 2))
    p = 2 * (1 - student_cdf(abs(t), n - 2))
    return rho, p


def student_cdf(t, df):
    x = df / (df + t * t)
    ib = betainc(df / 2.0, 0.5, x)
    return 1 - 0.5 * ib


def betainc(a, b, x):
    """Regularised incomplete beta, continued fraction."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) / a
    f, c, d = 1.0, 1.0, 0.0
    for i in range(0, 300):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        d = 1e-30 if abs(d) < 1e-30 else d
        d = 1.0 / d
        c = 1.0 + num / c
        c = 1e-30 if abs(c) < 1e-30 else c
        f *= c * d
        if abs(1 - c * d) < 1e-10:
            break
    result = front * (f - 1)
    if x > (a + 1) / (a + b + 2):
        return 1 - betainc(b, a, 1 - x)
    return result


def wilcoxon(d):
    """Paired signed-rank test, normal approximation with tie correction."""
    nz = [v for v in d if v != 0]
    n = len(nz)
    if n < 10:
        return float("nan"), float("nan"), n
    ranks = rankdata([abs(v) for v in nz])
    wp = sum(r for v, r in zip(nz, ranks) if v > 0)
    wm = sum(r for v, r in zip(nz, ranks) if v < 0)
    w = min(wp, wm)
    mean = n * (n + 1) / 4.0
    counts = {}
    for r in ranks:
        counts[r] = counts.get(r, 0) + 1
    tie = sum(c ** 3 - c for c in counts.values())
    var = (n * (n + 1) * (2 * n + 1) - tie / 2.0) / 24.0
    if var <= 0:
        return w, float("nan"), n
    z = (w - mean + 0.5) / math.sqrt(var)
    p = 2 * 0.5 * math.erfc(-(-abs(z)) / math.sqrt(2))
    return w, min(1.0, p), n


def report(label, rho, prho, w, pw, n, direction):
    print(f"{label}")
    print(f"  positions compared      : {n}")
    print(f"  Spearman rho            : {rho:.3f}   p = {fmt_p(prho)}")
    print(f"  Wilcoxon W              : {w:.1f}     p = {fmt_p(pw)}")
    print(f"  direction               : {direction}\n")


def fmt_p(p):
    if p != p:
        return "n/a"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.4f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--label-a", default="run A")
    ap.add_argument("--label-b", default="run B")
    ap.add_argument("--drop-flagged", action="store_true",
                    help="exclude low-confidence positions from either run")
    ap.add_argument("--by-domain", action="store_true")
    ap.add_argument("--csv")
    args = ap.parse_args()

    A, B = read_res(args.a), read_res(args.b)
    shared = sorted(set(A) & set(B))
    if not shared:
        sys.exit("no shared positions — are both files referenced to YFT67316.1?")
    if set(A) != set(B):
        print(f"! position sets differ: {len(A)} vs {len(B)}, "
              f"{len(shared)} shared\n", file=sys.stderr)

    for pos in shared:
        if A[pos]["aa"] != B[pos]["aa"]:
            sys.exit(f"residue mismatch at position {pos}: "
                     f"{A[pos]['aa']} vs {B[pos]['aa']} — different reference?")

    flagged = flagged_positions(A) | flagged_positions(B)
    used = [p for p in shared if not (args.drop_flagged and p in flagged)]

    print(f"\n{'=' * 70}")
    print("Cross-level comparison of site-specific evolutionary rates")
    print("=" * 70)
    print(f"A: {args.label_a}  ({args.a})")
    print(f"B: {args.label_b}  ({args.b})")
    print(f"Shared reference positions: {len(shared)}")
    print(f"Low-confidence positions  : {len(flagged)}"
          f"{' (excluded)' if args.drop_flagged else ' (retained)'}")
    print(f"Positions used            : {len(used)}\n")

    sa = [A[p]["score"] for p in used]
    sb = [B[p]["score"] for p in used]
    ra, rb = rankdata(sa), rankdata(sb)

    rho, prho = spearman(sa, sb)

    w_r, p_r, n_r = wilcoxon([x - y for x, y in zip(ra, rb)])
    shift = sum(1 for x, y in zip(ra, rb) if y > x)
    dir_r = (f"{shift} of {n_r} positions rank relatively faster in B, "
             f"{n_r - shift} slower")
    report("Ranks (primary — normalisation-invariant)",
           rho, prho, w_r, p_r, len(used), dir_r)

    w_s, p_s, n_s = wilcoxon([x - y for x, y in zip(sa, sb)])
    mean_d = sum(x - y for x, y in zip(sa, sb)) / len(used)
    report("Normalised scores (secondary — both runs are z-scored)",
           rho, prho, w_s, p_s, len(used),
           f"mean difference A-B = {mean_d:+.3f}")

    if args.by_domain:
        print("By domain (Spearman rho, A vs B):\n")
        print(f"{'Domain':<24}{'Residues':>10}{'n':>5}{'rho':>8}{'p':>12}")
        print("-" * 59)
        for name, start, end in DOMAINS:
            idx = [p for p in used if start <= p <= end]
            if len(idx) < 4:
                print(f"{name:<24}{f'{start}-{end}':>10}{len(idx):>5}"
                      f"{'n/a':>8}{'too few':>12}")
                continue
            r, pr = spearman([A[p]["score"] for p in idx],
                             [B[p]["score"] for p in idx])
            print(f"{name:<24}{f'{start}-{end}':>10}{len(idx):>5}"
                  f"{r:>8.3f}{fmt_p(pr):>12}")
        print()

    if args.csv:
        with open(args.csv, "w") as fh:
            fh.write("position,residue,score_A,score_B,rank_A,rank_B,"
                     "low_confidence\n")
            for i, p in enumerate(used):
                fh.write(f"{p},{A[p]['aa']},{sa[i]:.4f},{sb[i]:.4f},"
                         f"{ra[i]:.1f},{rb[i]:.1f},"
                         f"{'yes' if p in flagged else 'no'}\n")
        print(f"wrote {args.csv}\n")


if __name__ == "__main__":
    main()
