#!/usr/bin/env python3
"""
compare_sets_vs_3vop.py — compare predicted trimers against the vaccinia A27
crystal assembly (3VOP), one row per model.

Purpose: decide objectively how close each restraint set brings the model to
the crystal arrangement, and whether sets differ from one another. Written to
replace the interactive PyMOL comparison of S1, S2 and S4 with a saved table.

Usage, run from 7_structure_prediction/analysis/:

    # every model of every run in the screen
    python compare_sets_vs_3vop.py ../2_restraints_screening \\
        --xtal ../2_restraints_screening/3vop-assembly1.cif \\
        --out comparison_vs_3vop

    # named runs only
    python compare_sets_vs_3vop.py \\
        ../2_restraints_screening/results_MPXV_A29_S1-f_seed1 \\
        ../2_restraints_screening/results_MPXV_A29_S2-f_seed1 \\
        ../2_restraints_screening/results_MPXV_A29_S4-f_seed1 \\
        --xtal ../2_restraints_screening/3vop-assembly1.cif --out comparison_S1_S2_S4

The crystal is prepared as in the interactive session: non-polymer atoms and
zero-occupancy atoms removed, then renumbered by +20 so that its residues
carry native A27 numbering (21-84).

Only antiparallel models are comparable, so each model is first classified
from its helix axes and skipped if it is not antiparallel (that fact is
recorded in the table). Chains are then relabelled to the crystal convention
- A = parallel chain presenting the g face to the antiparallel chain,
B = parallel chain presenting the e face, C = the antiparallel chain - so that
chain A of a model is always compared with chain A of the crystal. Of the two
parallel chains, the one turning its g face towards the antiparallel chain
becomes A and the other B; the chain map applied to each model is recorded in
the table.

For every model and each residue window it reports:
  assembly    RMSD with all three chains fitted together
  chain_X     RMSD of chain X fitted on its own (fold accuracy, not packing)
  on_X:Y      deviation of chain Y with the structure fitted on chain X only
              (packing accuracy: how far the other chains sit once one is
              anchored)

Requires numpy and gemmi.
"""

import argparse
import csv
import itertools
import sys
from pathlib import Path

import numpy as np

try:
    import gemmi
except ImportError:
    sys.exit("gemmi is required:  pip install gemmi")

# ---------------------------------------------------------------- settings

WINDOWS = [(45, 84), (47, 65)]   # residue ranges compared, native numbering
AXIS_RANGE = (47, 65)            # residues used for the helix-axis fit
E_FACE = {48, 55, 62}            # heptad e positions, native numbering
G_FACE = {50, 57, 64}            # heptad g positions
FACE_CUTOFF = 8.0                # Å, Cβ-Cβ, for deciding which face points in
XTAL_SHIFT = 20                  # added to crystal numbering (1-64 -> 21-84)
PARALLEL_DOT = 0.7               # axis dot product above this: same direction


# ---------------------------------------------------------------- structures

def load_model(path, shift=0, drop_zero_occupancy=False):
    """Return {(chain, residue): {atom: coord}} for one structure."""
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_ligands_and_waters()
    st.setup_entities()
    out = {}
    for chain in st[0]:
        for res in chain:
            num = res.seqid.num + shift
            atoms = {}
            for at in res:
                if drop_zero_occupancy and at.occ < 0.01:
                    continue
                atoms[at.name] = np.array([at.pos.x, at.pos.y, at.pos.z])
            if atoms:
                out[(chain.name, num)] = atoms
    return out


def ca(struct, chain, lo, hi):
    """CA coordinates of one chain over a residue window, in residue order."""
    nums = sorted(n for (c, n) in struct if c == chain and lo <= n <= hi)
    return nums, np.array([struct[(chain, n)]["CA"] for n in nums
                           if "CA" in struct[(chain, n)]])


def cbeta(struct, chain, num):
    """Cβ of one residue, falling back to CA for glycine."""
    atoms = struct.get((chain, num))
    if not atoms:
        return None
    return atoms.get("CB", atoms.get("CA"))


def chains_of(struct):
    return sorted({c for (c, _) in struct})


# ---------------------------------------------------------------- geometry

def axis(struct, chain):
    """Unit vector along the helix, pointing from low to high residue number."""
    _, P = ca(struct, chain, *AXIS_RANGE)
    if len(P) < 8:
        return None
    v = np.linalg.svd(P - P.mean(axis=0), full_matrices=False)[2][0]
    if np.dot(P[-1] - P[0], v) < 0:
        v = -v
    return v / np.linalg.norm(v)


def classify(struct):
    """Return (class, antiparallel chain). Class is P, AP or X.

    AP means one chain runs opposite to the other two, as in the crystal.
    """
    ax = {c: axis(struct, c) for c in chains_of(struct)}
    if any(v is None for v in ax.values()) or len(ax) != 3:
        return "X", None
    dots = {f"{a}{b}": float(np.dot(ax[a], ax[b]))
            for a, b in itertools.combinations(sorted(ax), 2)}
    v = list(dots.values())
    if all(x > PARALLEL_DOT for x in v):
        return "P", None
    if sum(x > PARALLEL_DOT for x in v) == 1 and sum(x < -PARALLEL_DOT for x in v) == 2:
        anti = [c for c in ax
                if all(np.dot(ax[c], ax[o]) < 0 for o in ax if o != c)]
        if len(anti) == 1:
            return "AP", anti[0]
    return "X", None


def face_score(struct, parallel_chain, anti_chain):
    """How strongly the parallel chain turns its g face, rather than its e face,
    towards the antiparallel chain.

    Positive means g, negative means e. The count of face residues within the
    cutoff decides; where the counts tie, the mean closest approach of each
    face is used instead, so that a decision is still possible when the two
    chains sit further apart than the cutoff.
    """
    anti_cb = [cbeta(struct, anti_chain, n)
               for n in range(AXIS_RANGE[0], AXIS_RANGE[1] + 1)]
    anti_cb = np.array([p for p in anti_cb if p is not None])
    if len(anti_cb) == 0:
        return None
    dist = {}
    for r in sorted(E_FACE | G_FACE):
        p = cbeta(struct, parallel_chain, r)
        if p is not None:
            dist[r] = float(np.linalg.norm(anti_cb - p, axis=1).min())
    de = [dist[r] for r in E_FACE if r in dist]
    dg = [dist[r] for r in G_FACE if r in dist]
    if not de or not dg:
        return None
    ne = sum(1 for d in de if d < FACE_CUTOFF)
    ng = sum(1 for d in dg if d < FACE_CUTOFF)
    if ng != ne:
        return float(ng - ne)
    return float(np.mean(de) - np.mean(dg)) / 100.0


def relabel(struct, anti_chain):
    """Rename chains to the crystal convention: A = g face, B = e face, C = anti."""
    others = [c for c in chains_of(struct) if c != anti_chain]
    scores = {c: face_score(struct, c, anti_chain) for c in others}
    if any(v is None for v in scores.values()) or len(others) != 2:
        return None, scores
    g_chain = max(others, key=lambda c: scores[c])
    e_chain = [c for c in others if c != g_chain][0]
    if scores[g_chain] == scores[e_chain]:
        return None, scores
    mapping = {g_chain: "A", e_chain: "B", anti_chain: "C"}
    return ({(mapping[c], n): atoms for (c, n), atoms in struct.items()},
            mapping)


def kabsch_rmsd(P, Q):
    """RMSD after optimal superposition of P onto Q."""
    Pc, Qc = P - P.mean(axis=0), Q - Q.mean(axis=0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1.0, 1.0, d]) @ Vt
    return float(np.sqrt(((Pc @ R - Qc) ** 2).sum(axis=1).mean()))


def kabsch_transform(P, Q):
    """Rotation and translation putting P onto Q."""
    Pm, Qm = P.mean(axis=0), Q.mean(axis=0)
    U, _, Vt = np.linalg.svd((P - Pm).T @ (Q - Qm))
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1.0, 1.0, d]) @ Vt
    return R, Qm - Pm @ R


def paired(model, xtal, chain, lo, hi):
    """CA coordinates present in both structures, in the same residue order."""
    nums = [n for n in range(lo, hi + 1)
            if (chain, n) in model and (chain, n) in xtal
            and "CA" in model[(chain, n)] and "CA" in xtal[(chain, n)]]
    return (np.array([model[(chain, n)]["CA"] for n in nums]),
            np.array([xtal[(chain, n)]["CA"] for n in nums]))


# ---------------------------------------------------------------- comparison

def compare(model, xtal, lo, hi):
    """All measurements for one model over one residue window."""
    pairs = {c: paired(model, xtal, c, lo, hi) for c in "ABC"}
    if any(len(p[0]) < 5 for p in pairs.values()):
        return None

    row = {}
    P = np.vstack([pairs[c][0] for c in "ABC"])
    Q = np.vstack([pairs[c][1] for c in "ABC"])
    row["assembly"] = kabsch_rmsd(P, Q)
    row["n_ca"] = len(P)

    for c in "ABC":
        row[f"chain_{c}"] = kabsch_rmsd(*pairs[c])

    for anchor in "ABC":
        R, t = kabsch_transform(pairs[anchor][0], pairs[anchor][1])
        for c in "ABC":
            moved = pairs[c][0] @ R + t
            dev = float(np.sqrt(((moved - pairs[c][1]) ** 2).sum(axis=1).mean()))
            row[f"on_{anchor}:{c}"] = dev
        row[f"on_{anchor}_others"] = max(
            row[f"on_{anchor}:{c}"] for c in "ABC" if c != anchor)
    return row


def set_name(path: Path) -> str:
    """Restraint set taken from the run folder name, e.g. ..._S2-f_seed1 -> S2-f."""
    for part in reversed(path.parts):
        if part.startswith("results_"):
            bits = part.split("_")
            return "_".join(b for b in bits if b.startswith("S") or b.startswith("seed"))
    return path.parent.name


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(
        description="Compare predicted trimers against the 3VOP crystal assembly.")
    ap.add_argument("inputs", nargs="+",
                    help="run folders, or a folder containing results_*/ folders")
    ap.add_argument("--xtal", required=True,
                    help="path to 3vop-assembly1.cif")
    ap.add_argument("--out", default="comparison_vs_3vop",
                    help="output folder (default: comparison_vs_3vop)")
    ap.add_argument("--shift", type=int, default=XTAL_SHIFT,
                    help="value added to crystal residue numbers (default 20)")
    args = ap.parse_args()

    xtal = load_model(args.xtal, shift=args.shift, drop_zero_occupancy=True)
    cls, anti = classify(xtal)
    if cls != "AP":
        sys.exit(f"The crystal did not classify as antiparallel ({cls}); check --shift.")
    xtal, mapping = relabel(xtal, anti)
    if xtal is None:
        sys.exit("Could not assign heptad faces in the crystal.")
    print(f"crystal: antiparallel chain {anti}, chains mapped {mapping}")

    cifs = []
    for item in args.inputs:
        p = Path(item)
        if p.is_file():
            cifs.append(p)
        else:
            cifs += sorted(p.glob("**/predictions/**/*model_*.cif"))
    cifs = sorted(set(cifs))
    if not cifs:
        sys.exit("No model files found.")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for cif in cifs:
        model = load_model(cif)
        cls, anti = classify(model)
        tag = set_name(cif)
        base = {"set": tag, "model": cif.name, "class": cls}
        if cls != "AP":
            rows.append({**base, "window": "-", "note": "not antiparallel, not compared"})
            print(f"  {tag:28s} {cif.name:34s} {cls}  skipped")
            continue
        relabelled, mapping = relabel(model, anti)
        if relabelled is None:
            rows.append({**base, "window": "-", "note": "heptad faces not assignable"})
            print(f"  {tag:28s} {cif.name:34s} AP  faces not assignable")
            continue
        for lo, hi in WINDOWS:
            res = compare(relabelled, xtal, lo, hi)
            if res is None:
                continue
            rows.append({**base, "window": f"{lo}-{hi}",
                         "chain_map": "".join(f"{k}->{v}" for k, v in sorted(mapping.items())),
                         **{k: (round(v, 2) if isinstance(v, float) else v)
                            for k, v in res.items()}})
        a = next(r for r in rows[::-1] if r.get("window") == f"{WINDOWS[0][0]}-{WINDOWS[0][1]}")
        print(f"  {tag:28s} {cif.name:34s} AP  assembly {a['assembly']:.2f} Å")

    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    table = out / "per_model.tsv"
    with table.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)

    # summary: range over the models of each set
    summary = out / "per_set.tsv"
    keys = ["assembly", "on_A_others", "on_B_others", "on_C_others"]
    with summary.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["set", "window", "n_models_AP"]
                   + [f"{k}_min" for k in keys] + [f"{k}_max" for k in keys])
        groups = {}
        for r in rows:
            if r.get("window") in (None, "-"):
                continue
            groups.setdefault((r["set"], r["window"]), []).append(r)
        for (tag, window), rs in sorted(groups.items()):
            vals = [[r[k] for r in rs if k in r] for k in keys]
            w.writerow([tag, window, len(rs)]
                       + [f"{min(v):.2f}" if v else "" for v in vals]
                       + [f"{max(v):.2f}" if v else "" for v in vals])

    print(f"\nWrote {table} and {summary}")


if __name__ == "__main__":
    main()
