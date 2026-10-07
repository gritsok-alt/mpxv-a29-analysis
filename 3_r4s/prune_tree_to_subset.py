#!/usr/bin/env python3
"""Prune the 60-taxon A29L ML tree down to the taxa present in a subset alignment.

Produces the fixed topology for the second Rate4Site run (Orthopoxvirus level).
Branch lengths are preserved: when a tip is removed, Bio.Phylo collapses the
resulting degree-two node and adds the two branch lengths together, so
root-to-tip distances among the retained taxa are unchanged.

The topology is pruned from the 60-taxon maximum-likelihood tree and is NOT
re-inferred, which is what the Methods state and what keeps the two Rate4Site
levels comparable.

Usage:
    python3 prune_tree_to_subset.py \
        --tree A29L.treefile \
        --msa  14sequences_genus.fasta \
        --out  A29L_orthopox14.treefile
"""

import argparse
import io
import re
import sys
from pathlib import Path


def norm_id(x):
    """Strip a Jalview-style /start-end suffix, matching run_rate4site.py."""
    return re.sub(r"/\d+-\d+$", "", x)


def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(norm_id(line[1:].strip().split()[0]))
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True, help="60-taxon ML tree (Newick)")
    ap.add_argument("--msa", required=True, help="subset alignment (aligned FASTA)")
    ap.add_argument("--out", required=True, help="pruned tree to write")
    args = ap.parse_args()

    from Bio import Phylo

    keep = read_fasta_ids(args.msa)
    if not keep:
        sys.exit(f"no sequences found in {args.msa}")
    if len(set(keep)) != len(keep):
        sys.exit("duplicate identifiers in the subset alignment")

    tree = Phylo.read(args.tree, "newick")
    tips = {norm_id(t.name): t for t in tree.get_terminals() if t.name}

    missing = [k for k in keep if k not in tips]
    if missing:
        sys.exit("subset taxa absent from the tree:\n  " + "\n  ".join(missing))

    drop = [name for name in tips if name not in set(keep)]
    print(f"tree tips     : {len(tips)}")
    print(f"subset taxa   : {len(keep)}")
    print(f"pruning away  : {len(drop)}")

    for name in drop:
        tree.prune(tips[name])

    remaining = sorted(norm_id(t.name) for t in tree.get_terminals() if t.name)
    if remaining != sorted(keep):
        sys.exit("pruned tip set does not match the alignment — aborting")

    # Drop internal labels and support values; run_rate4site.py does this too,
    # but a clean tree here keeps the deposited intermediate readable.
    for cl in tree.get_nonterminals():
        cl.name = None
        cl.confidence = None

    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    nwk = buf.getvalue().strip()
    nwk = re.sub(r"(?<=:)([\d.]+[eE][-+]?\d+)",
                 lambda m: f"{float(m.group(1)):.10f}", nwk)

    with open(args.out, "w", newline="\n", encoding="ascii") as fh:
        fh.write(nwk + "\n")

    print(f"\nwrote {Path(args.out).name} with {len(remaining)} tips")
    print("tips retained:")
    for name in remaining:
        print(f"  {name}")


if __name__ == "__main__":
    main()
