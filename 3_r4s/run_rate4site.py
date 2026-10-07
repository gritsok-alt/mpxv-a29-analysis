#!/usr/bin/env python3
"""
Rate4Site run for MPXV A29 (OPG154), referenced to YFT67316.1.

Reproduces the primary run documented in README_rate4site.md:

    rate4site -s aln.fasta -t tree.nwk -a T60 -Mj -ib -k 16 \
              -o r60_jtt.res -y r60_jtt.unnorm -x r60_jtt.tree

Two preprocessing steps are necessary and are handled here, because Rate4Site's
Newick parser fails on both silently:

  1. Taxon labels containing '|' are renamed to T01..Tnn. A pipe inside a name
     is passed to a shell before Rate4Site sees it, so a label such as
     'Orthopoxvirus|YFT67316.1|...' is parsed as three piped commands.
  2. IQ-TREE writes internal-node labels of the form '75.1/0.856/94', which
     Rate4Site cannot parse. These are stripped; branch lengths and topology
     are preserved.

A name map is written so the short names can be restored afterwards.

Usage:
    python run_rate4site.py --msa aln.fasta --tree A29L.treefile --ref YFT67316.1
"""

import argparse
import io
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


def norm_id(x):
    """Strip a Jalview-style /start-end suffix."""
    return re.sub(r"/\d+-\d+$", "", x)


def sanitise_alignment(msa_path, out_path):
    """Rename taxa to T01..Tnn; return {original: short} and the record order."""
    from Bio import AlignIO

    aln = AlignIO.read(msa_path, "fasta")
    for r in aln:
        r.id = norm_id(r.id)
        r.description = ""
    name_map = {r.id: f"T{i:02d}" for i, r in enumerate(aln, 1)}
    for r in aln:
        r.id = name_map[norm_id(r.id)]
        r.name = r.id
    AlignIO.write(aln, out_path, "fasta")
    print(f"  alignment: {len(aln)} sequences x {aln.get_alignment_length()} columns")
    return name_map


def sanitise_tree(tree_path, out_path, name_map):
    """Strip internal labels, rename tips, write LF-terminated plain-decimal Newick."""
    from Bio import Phylo

    tree = Phylo.read(tree_path, "newick")
    for cl in tree.get_nonterminals():
        cl.name = None
        cl.confidence = None
    missing = []
    for tip in tree.get_terminals():
        key = norm_id(tip.name)
        if key not in name_map:
            missing.append(key)
        tip.name = name_map.get(key, key)
    if missing:
        sys.exit(f"tree tips absent from the alignment: {missing[:5]}")

    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    nwk = buf.getvalue().strip()
    # Rate4Site cannot read scientific notation in branch lengths
    nwk = re.sub(r"(?<=:)([\d.]+[eE][-+]?\d+)",
                 lambda m: f"{float(m.group(1)):.10f}", nwk)
    with open(out_path, "w", newline="\n", encoding="ascii") as fh:
        fh.write(nwk + "\n")

    residual = re.findall(r"\)([^:,();]+)", nwk)
    print(f"  tree: {len(tree.get_terminals())} tips, "
          f"residual internal labels: {residual[:3] or 'none'}")
    return tree


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--msa", required=True, help="aligned FASTA")
    ap.add_argument("--tree", required=True, help="ML tree (Newick)")
    ap.add_argument("--ref", default="YFT67316.1",
                    help="accession of the reference sequence (default %(default)s)")
    ap.add_argument("--out", default="r4s_run", help="output directory")
    ap.add_argument("--tag", default="r60_jtt", help="output basename")
    ap.add_argument("--matrix", default="-Mj",
                    help="substitution matrix flag: -Mj JTT, -Mw WAG, -Ml LG")
    ap.add_argument("--categories", type=int, default=16,
                    help="discrete gamma categories (default 16)")
    ap.add_argument("--no-bl-opt", action="store_true",
                    help="do not re-optimise branch lengths (-bn)")
    ap.add_argument("--exe", default="rate4site")
    args = ap.parse_args()

    if not shutil.which(args.exe):
        sys.exit(f"{args.exe} not found on PATH.\n"
                 "Build it from https://github.com/barakav/r4s_for_collab "
                 "or install via conda: conda install -c bioconda rate4site")

    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    print("Preprocessing")
    aln_c = out / "aln.fasta"
    tree_c = out / "tree.nwk"
    name_map = sanitise_alignment(args.msa, aln_c)
    sanitise_tree(args.tree, tree_c, name_map)

    hits = [orig for orig in name_map if args.ref in orig]
    if len(hits) != 1:
        sys.exit(f"reference {args.ref} matched {len(hits)} sequences: {hits}")
    ref_short = name_map[hits[0]]
    print(f"  reference: {hits[0]}  ->  {ref_short}")

    (out / "namemap.json").write_text(json.dumps(name_map, indent=2))

    cmd = [args.exe,
           "-s", str(aln_c),
           "-t", str(tree_c),
           "-a", ref_short,
           args.matrix, "-ib",
           "-k", str(args.categories),
           "-o", str(out / f"{args.tag}.res"),
           "-y", str(out / f"{args.tag}.unnorm"),
           "-x", str(out / f"{args.tag}.tree")]
    if args.no_bl_opt:
        cmd.append("-bn")

    print("\nRunning\n  " + " ".join(cmd))
    r = subprocess.run(cmd, cwd=out, capture_output=True, text=True)
    (out / f"{args.tag}.log").write_text((r.stdout or "") + "\n" + (r.stderr or ""))

    res = out / f"{args.tag}.res"
    if not (res.exists() and res.stat().st_size):
        print((r.stdout or "")[-1200:], (r.stderr or "")[-600:])
        sys.exit(f"\nfailed - see {args.tag}.log")

    n = sum(1 for l in res.read_text().splitlines()
            if re.match(r"\s*\d+\s+\w\s+-?[\d.]+", l))
    alpha = next((l.split()[-1] for l in res.read_text().splitlines()
                  if "alpha parameter" in l), "?")
    print(f"\n  {res.name}: {n} positions, fitted gamma alpha = {alpha}")
    print(f"\nWritten to {out}/")
    print(f"  {args.tag}.res      rates, in reference coordinates")
    print(f"  {args.tag}.unnorm   un-normalised rates")
    print(f"  {args.tag}.tree     tree with re-optimised branch lengths")
    print(f"  aln.fasta, tree.nwk, namemap.json, {args.tag}.log")
    print(f"\nNext:  python plot_rate4site_A29.py --res {out}/{args.tag}.res")


if __name__ == "__main__":
    main()
