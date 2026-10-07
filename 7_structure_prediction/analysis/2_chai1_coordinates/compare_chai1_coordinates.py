#!/usr/bin/env python3
"""
compare_models.py — compare the coordinates of predicted structures.

Purpose: decide objectively whether two or more predictions are the same
structure (identical coordinates), a rigid-body copy of each other, or
genuinely different models. Written to document why the Chai-1 runs were
excluded for returning identical coordinates.

Usage, from the folder holding the runs to compare:

    python compare_models.py                 # compares every .cif/.pdb below the current folder
    python compare_models.py <folder> ...    # or compare specific folders/files
    python compare_models.py --rank 0        # only files whose name contains "model_0"/"rank_0"

Requires numpy and biopython only.

For every pair of structures it reports:
  max_dev   the largest distance between corresponding atoms, as the files
            are written (no superposition). 0.000 means byte-identical
            coordinates; a small value means the same model rewritten.
  rmsd_sup  RMSD after optimal superposition (Kabsch). If max_dev is large
            but rmsd_sup is ~0, the models are the same structure in a
            different frame of reference.
Pairs are then grouped: structures whose rmsd_sup is below the tolerance are
reported as one group of duplicates.
"""

import sys
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter("ignore", PDBConstructionWarning)

TOL = 0.10          # Å; below this, two structures are called identical
ATOM = "CA"         # atom used for the comparison


def load(path: Path):
    """Return {(chain, residue_number): coordinate} for the chosen atom."""
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() == ".cif" else PDBParser(QUIET=True)
    model = next(parser.get_structure(path.stem, str(path)).get_models())
    coords = {}
    for chain in model:
        for res in chain:
            if ATOM in res:
                coords[(chain.id, res.id[1])] = res[ATOM].get_coord()
    return coords


def label(path: Path) -> str:
    """Short identifier: parent folder + file name, so that runs sharing a
    file name (pred.rank_0.cif) remain distinguishable."""
    return f"{path.parent.name}/{path.name}"


def common(a, b):
    keys = sorted(set(a) & set(b))
    return (np.array([a[k] for k in keys]),
            np.array([b[k] for k in keys]),
            len(keys))


def kabsch_rmsd(P, Q):
    """RMSD after optimal superposition of P onto Q."""
    P = P - P.mean(axis=0)
    Q = Q - Q.mean(axis=0)
    V, S, W = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ W))
    D = np.diag([1.0, 1.0, d])
    P = P @ (V @ D @ W)
    return float(np.sqrt(((P - Q) ** 2).sum() / len(P)))


def main(argv):
    rank = None
    args = []
    for i, a in enumerate(argv):
        if a == "--rank":
            rank = argv[i + 1]
        elif argv[i - 1] != "--rank":
            args.append(a)

    roots = [Path(a) for a in args] or [Path(".")]
    files = []
    for root in roots:
        if root.is_file():
            files.append(root)
        else:
            files += [p for p in root.rglob("*") if p.suffix.lower() in (".cif", ".pdb")]
    if rank is not None:
        files = [f for f in files if f"model_{rank}" in f.name or f"rank_{rank}" in f.name]
    files = sorted(set(files))

    if len(files) < 2:
        sys.exit(f"Need at least two structures to compare; found {len(files)}.")

    print(f"{len(files)} structures, comparing {ATOM} atoms\n")
    structures = {}
    for f in files:
        c = load(f)
        structures[f] = c
        print(f"  {len(c):4d} residues  {f}")

    print(f"\n{'file A':<45}{'file B':<45}{'n':>5}{'max_dev':>10}{'rmsd_sup':>10}  verdict")
    print("-" * 125)
    identical = []
    for a, b in combinations(files, 2):
        P, Q, n = common(structures[a], structures[b])
        if n == 0:
            print(f"{label(a):<45}{label(b):<45}{0:>5}{'-':>10}{'-':>10}  no shared residues")
            continue
        max_dev = float(np.linalg.norm(P - Q, axis=1).max())
        rmsd = kabsch_rmsd(P, Q)
        if max_dev < 1e-3:
            verdict = "IDENTICAL coordinates"
        elif rmsd < TOL:
            verdict = "same structure, different frame"
        else:
            verdict = "different"
        if rmsd < TOL:
            identical.append((a, b))
        print(f"{label(a):<45}{label(b):<45}{n:>5}{max_dev:>10.3f}{rmsd:>10.3f}  {verdict}")

    # group duplicates
    groups = []
    for a, b in identical:
        for g in groups:
            if a in g or b in g:
                g.update({a, b})
                break
        else:
            groups.append({a, b})

    if groups:
        print("\nGroups of structures that are the same within "
              f"{TOL} Å RMSD after superposition:")
        for i, g in enumerate(groups, 1):
            print(f"  group {i}:")
            for f in sorted(g):
                print(f"    {f}")
    else:
        print(f"\nNo two structures agree within {TOL} Å RMSD.")


if __name__ == "__main__":
    main(sys.argv[1:])
