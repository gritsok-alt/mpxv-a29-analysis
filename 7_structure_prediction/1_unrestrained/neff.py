#!/usr/bin/env python3
"""
Effective alignment depth (Neff) of one or more MSAs.

Neff is the sum of sequence weights, each sequence weighted by the reciprocal
of the number of sequences in the alignment that resemble it at or above the
identity threshold. A cluster of near-identical sequences therefore contributes
approximately one effective sequence rather than many.

The threshold matters: a higher value clusters less aggressively and gives a
larger Neff, so figures computed at different thresholds are not comparable.
Report the threshold alongside every value.

Usage:
    python neff.py alignment.a3m [more.a3m ...] [-t 0.62] [-t 0.80]
    python neff.py <folder>                      # every .a3m below it

Handles a3m (lower-case insertions removed), aligned FASTA and Stockholm.
"""

import argparse
import sys
from pathlib import Path

import numpy as np


def read_alignment(path):
    """Return a list of aligned sequences, insertions stripped for a3m."""
    text = Path(path).read_text()

    if text.lstrip().startswith("# STOCKHOLM"):
        seqs = {}
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith(("#", "//")):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                seqs[parts[0]] = seqs.get(parts[0], "") + parts[1]
        raw = list(seqs.values())
    else:
        raw, cur = [], []
        for line in text.splitlines():
            if line.startswith(">"):
                if cur:
                    raw.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
        if cur:
            raw.append("".join(cur))

    # a3m: lower case marks an insertion relative to the query; drop those
    cleaned = ["".join(c for c in s if not c.islower() and c != ".") for s in raw]
    cleaned = [s for s in cleaned if s]
    if not cleaned:
        return []

    L = len(cleaned[0])
    keep = [s for s in cleaned if len(s) == L]
    if len(keep) < len(cleaned):
        print(f"    note: {len(cleaned) - len(keep)} sequence(s) of unequal "
              f"length discarded", file=sys.stderr)
    return keep


def neff(seqs, threshold, chunk=2000):
    """Sum of 1/n_i, where n_i counts sequences within `threshold` identity."""
    if not seqs:
        return None, 0, 0

    arr = np.frombuffer("".join(seqs).encode(), dtype="S1").reshape(len(seqs), -1)

    # restrict to positions where the query has a residue
    query_has_residue = arr[0] != b"-"
    A = arr[:, query_has_residue]
    n, L = A.shape
    if L == 0:
        return None, n, 0

    # integer encoding is faster to compare
    codes = np.unique(A)
    lookup = {c: i for i, c in enumerate(codes)}
    Ai = np.vectorize(lookup.get)(A).astype(np.int16)

    counts = np.zeros(n, dtype=np.int32)
    for start in range(0, n, chunk):
        block = Ai[start:start + chunk]
        # fraction of identical positions against every sequence
        ident = (block[:, None, :] == Ai[None, :, :]).sum(axis=2) / L
        counts[start:start + chunk] = (ident >= threshold).sum(axis=1)

    weights = 1.0 / np.maximum(counts, 1)
    return float(weights.sum()), n, L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("-t", "--threshold", type=float, action="append",
                    help="identity threshold, repeatable (default 0.62 and 0.80)")
    args = ap.parse_args()

    thresholds = args.threshold or [0.62, 0.80]

    files = []
    for item in args.inputs:
        p = Path(item).expanduser()
        if p.is_dir():
            files += sorted(p.rglob("*.a3m")) + sorted(p.rglob("*.sto")) \
                   + sorted(p.rglob("*.fasta")) + sorted(p.rglob("*.fa"))
        elif p.is_file():
            files.append(p)
    if not files:
        sys.exit("No alignment files found.")

    header = f"{'alignment':<38}{'seqs':>7}{'L':>6}"
    for t in thresholds:
        header += f"{f'Neff@{t:.0%}':>12}{f'Neff/L':>9}"
    print(header)
    print("-" * len(header))

    for f in files:
        seqs = read_alignment(f)
        if not seqs:
            print(f"{f.name:<38}   could not parse")
            continue
        line = f"{f.name[:37]:<38}"
        first = True
        for t in thresholds:
            v, n, L = neff(seqs, t)
            if first:
                line += f"{n:>7}{L:>6}"
                first = False
            line += f"{v:>12.2f}{v / L:>9.3f}" if v else f"{'-':>12}{'-':>9}"
        print(line)

    print(f"\nThresholds compared: {', '.join(f'{t:.0%}' for t in thresholds)}")
    print("Values at different thresholds are not comparable; quote the "
          "threshold with every figure.")


if __name__ == "__main__":
    main()
