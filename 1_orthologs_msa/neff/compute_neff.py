#!/usr/bin/env python3
"""
compute_neff.py — read every aln_nr*.fasta in this folder and print one row each.

Run from Project_2/N_eff/ :      python compute_neff.py

Needs numpy + biopython only. No MAFFT, no realignment, no arguments.
Any alignment file that is missing is simply skipped, so you can run this
before or after producing aln_nr99.fasta and aln_nr100.fasta.
"""

import csv
from pathlib import Path

import numpy as np
from Bio import SeqIO

THR  = 0.62               # two sequences are "neighbours" at >= 62 % identity
REF  = "YFT67316.1"       # MPXV clade I reference, defines core numbering
CORE = (44, 110)          # confidently aligned core, reference coordinates
L    = 110                # reference length, for N_eff / L

ORDER = ["aln_nr90.fasta", "aln_nr95.fasta", "aln_nr98.fasta",
         "aln_nr99.fasta", "aln_nr100.fasta"]


def read(f):
    d = {r.id: str(r.seq).upper() for r in SeqIO.parse(f, "fasta")}
    widths = {len(v) for v in d.values()}
    if len(widths) != 1:
        raise SystemExit(f"{f.name} is not aligned — column counts {sorted(widths)}")
    return d


def core_cols(aln, lo, hi):
    key = next((k for k in aln if REF in k), None)
    if key is None:
        raise SystemExit(f"reference {REF} not found in the alignment")
    cols, pos = [], 0
    for c, ch in enumerate(aln[key]):
        if ch not in "-.":
            pos += 1
            if lo <= pos <= hi:
                cols.append(c)
    return cols


def components(adj):
    """Connected components of the neighbour graph. N_eff can never be smaller
    than this, and equals it only when every component is fully connected."""
    n, seen, k = len(adj), set(), 0
    for s in range(n):
        if s in seen:
            continue
        k += 1
        stack = [s]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            stack += [v for v in np.flatnonzero(adj[u]) if v not in seen]
    return k


def stats(aln, cols=None):
    """N_eff = sum_i 1 / (number of sequences within 62 % identity of i,
    counting i itself). Identity is measured only over columns where BOTH
    sequences have a residue."""
    M = np.array([list(v) for v in aln.values()])
    if cols is not None:
        M = M[:, cols]
    gap = np.isin(M, list("-.")) | (M == "X")

    n = len(M)
    P = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            both = (~gap[i]) & (~gap[j])
            P[i, j] = P[j, i] = (M[i][both] == M[j][both]).mean() if both.sum() else 0.0

    adj = P >= THR
    w = 1.0 / adj.sum(1)
    iu = np.triu_indices(n, 1)
    used = ~np.all(gap, axis=0)            # columns that are not all-gap
    return dict(n=n, cols=int(used.sum()), neff=float(w.sum()),
                per_L=float(w.sum() / L), clusters=components(adj),
                edges=int(adj[iu].sum()), singletons=int((adj.sum(1) == 1).sum()),
                mean_id=float(P[iu].mean()), median_id=float(np.median(P[iu])))


here = Path(__file__).parent
rows = []

for name in ORDER:
    f = here / name
    if not f.exists():
        print(f"skipping {name} (not found)")
        continue
    rows.append((name.replace("aln_", "").replace(".fasta", ""), stats(read(f))))

f98 = here / "aln_nr98.fasta"
if f98.exists():
    a = read(f98)
    rows.append((f"nr98 core {CORE[0]}-{CORE[1]}", stats(a, cols=core_cols(a, *CORE))))

print(f"\n{'set':<20}{'n':>5}{'cols':>6}{'N_eff':>9}{'N_eff/L':>9}"
      f"{'clusters':>10}{'edges':>7}{'alone':>7}{'meanID':>8}{'medID':>8}")
print("-" * 89)
for name, s in rows:
    print(f"{name:<20}{s['n']:>5}{s['cols']:>6}{s['neff']:>9.2f}{s['per_L']:>9.3f}"
          f"{s['clusters']:>10}{s['edges']:>7}{s['singletons']:>7}"
          f"{s['mean_id']:>8.3f}{s['median_id']:>8.3f}")

with open(here / "neff_table.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["set", "n_sequences", "columns", "N_eff", "N_eff_per_L",
                "clusters", "edges", "singletons", "mean_identity", "median_identity"])
    for name, s in rows:
        w.writerow([name, s["n"], s["cols"], round(s["neff"], 3), round(s["per_L"], 4),
                    s["clusters"], s["edges"], s["singletons"],
                    round(s["mean_id"], 4), round(s["median_id"], 4)])

print(f"\nwrote neff_table.csv")
print("""
CHECK THESE FOUR NUMBERS against the manuscript — they must not have moved:
    nr90 18.01     nr95 17.91     nr98 17.85     core 13.89
If they have, something went wrong; stop and say so.

WHAT THE COLUMNS MEAN
  N_eff      effective number of independent sequences
  clusters   groups of sequences linked by >= 62 % identity. N_eff >= clusters.
  edges      how many sequence pairs exceed 62 % identity
  alone      sequences with no neighbour at all; each contributes exactly 1.0
""")
