#!/usr/bin/env python3
"""
Relative solvent accessibility across the restrained-model ensemble.

Single-model accessibility is unreliable for residues in the C-terminal half
of this protein: the ensemble RMSF rises to 3-4 A by residue 100, so a chain
whose terminus happens to fold across the bundle in one model will bury
residues that are exposed in another. This script therefore reports the
distribution across the accepted ensemble rather than one model's value.

Acceptance matches the ensemble analysis:
  - antiparallel topology (the one all-parallel model is excluded)
  - ipTM >= MIN_IPTM       (two further models fall below the threshold)

Accessibility is computed with FreeSASA if available, otherwise with
Shrake-Rupley via Biopython. Relative values use the theoretical maxima of
Tien et al. (2013) PLoS ONE 8:e80635.

Usage:
    python ensemble_sasa.py <results_dir> [--out <dir>] [--min-iptm 0.60]
"""

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

try:
    import gemmi
except ImportError:
    sys.exit("gemmi is required:  pip install gemmi")


# ---------------------------------------------------------------- configuration

FIRST_RESIDUE = 1
CC_START, CC_END = 44, 100

# Residue groups to report, with what the literature says about each
GROUPS = [
    ("Hydrophobic core (a/d positions)", [47, 51, 54, 58, 61, 65],
     "required for self-assembly; expected buried"),
    ("A18-binding surface", [81, 87, 88, 94, 95, 98, 99],
     "alanine substitution abolishes A17/A18 binding; expected exposed"),
    ("Cysteine pair", [71, 72],
     "intermolecular disulfide to A28; Cys71 at g, Cys72 at a"),
    ("Polar core layer", [75],
     "Asn75 at d; the only buried polar position"),
]

# Tien et al. (2013), theoretical maximum ASA per residue (A^2)
MAX_ASA = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
    "GLU": 223, "GLN": 225, "GLY": 104, "HIS": 224, "ILE": 197,
    "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
    "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174,
}

BURIED_BELOW = 20.0     # % relative accessibility
EXPOSED_ABOVE = 25.0


# ---------------------------------------------------------------- SASA backend

def sasa_freesasa(path):
    """Per-residue absolute SASA via FreeSASA. Returns {(chain, resnum): area}."""
    import freesasa
    s = freesasa.Structure(str(path))
    r = freesasa.calc(s)
    out = {}
    for i in range(s.nAtoms()):
        key = (s.chainLabel(i).strip(), int(s.residueNumber(i)))
        out[key] = out.get(key, 0.0) + r.atomArea(i)
    return out


def sasa_biopython(path):
    from Bio.PDB import MMCIFParser, PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    parser = (MMCIFParser(QUIET=True) if str(path).endswith(".cif")
              else PDBParser(QUIET=True))
    st = parser.get_structure("m", str(path))
    ShrakeRupley().compute(st[0], level="R")
    out = {}
    for ch in st[0]:
        for res in ch:
            out[(ch.id.strip(), res.id[1])] = res.sasa
    return out


def pick_backend():
    try:
        import freesasa  # noqa: F401
        return sasa_freesasa, "FreeSASA"
    except ImportError:
        pass
    try:
        import Bio  # noqa: F401
        return sasa_biopython, "Biopython Shrake-Rupley"
    except ImportError:
        sys.exit("Install one of:\n  pip install freesasa\n  pip install biopython")


# ---------------------------------------------------------------- model handling

def helix_axis(chain):
    pts = []
    for r in chain:
        n = r.seqid.num + FIRST_RESIDUE - 1
        if CC_START <= n <= CC_END:
            a = r.find_atom("CA", "*")
            if a:
                pts.append([a.pos.x, a.pos.y, a.pos.z])
    if len(pts) < 8:
        return None
    P = np.array(pts)
    v = np.linalg.svd(P - P.mean(0), full_matrices=False)[2][0]
    if np.dot(P[-1] - P[0], v) < 0:
        v = -v
    return v / np.linalg.norm(v)


def classify_and_map(st):
    """
    Returns (class, {old_chain: new_label}) with C the antiparallel chain,
    A the chain presenting the heptad g face to C, B the other.
    """
    axes = {ch.name: helix_axis(ch) for ch in st[0]}
    axes = {k: v for k, v in axes.items() if v is not None}
    if len(axes) < 3:
        return "?", None

    dots = [np.dot(axes[a], axes[b])
            for a, b in itertools.combinations(sorted(axes), 2)]
    if all(d > 0.7 for d in dots):
        return "P", None
    if not (sum(d > 0.7 for d in dots) == 1 and sum(d < -0.7 for d in dots) == 2):
        return "X", None

    names = sorted(axes)
    anti = [n for n in names
            if all(np.dot(axes[n], axes[m]) < 0 for m in names if m != n)]
    if len(anti) != 1:
        return "X", None
    c = anti[0]
    others = [n for n in names if n != c]

    E, G = {48, 55, 62}, {50, 57, 64}

    def cb(ch, nums):
        pts = []
        for r in st[0][ch]:
            if r.seqid.num + FIRST_RESIDUE - 1 in nums:
                a = r.find_atom("CB", "*") or r.find_atom("CA", "*")
                if a:
                    pts.append(([a.pos.x, a.pos.y, a.pos.z],
                                r.seqid.num + FIRST_RESIDUE - 1))
        return pts

    def face(par):
        Q = np.array([[a.pos.x, a.pos.y, a.pos.z] for r in st[0][c]
                      for a in [r.find_atom("CB", "*") or r.find_atom("CA", "*")]
                      if a and CC_START <= r.seqid.num + FIRST_RESIDUE - 1 <= CC_END])
        hits = set()
        for xyz, num in cb(par, E | G):
            if len(Q) and np.min(np.linalg.norm(Q - np.array(xyz), axis=1)) <= 8.0:
                hits.add(num)
        return len(hits & E), len(hits & G)

    scored = [(par, face(par)) for par in others]
    g_chain = max(scored, key=lambda t: t[1][1] - t[1][0])[0]
    e_chain = [n for n in others if n != g_chain][0]
    return "AP", {g_chain: "A", e_chain: "B", c: "C"}


def relabel_and_write(st, mapping, dest, keep=None):
    tmp = {old: f"~{new}" for old, new in mapping.items()}
    for ch in st[0]:
        if ch.name in tmp:
            ch.name = tmp[ch.name]
    for ch in st[0]:
        if ch.name.startswith("~"):
            ch.name = ch.name[1:]

    if keep is not None:
        lo, hi = keep
        for ch in st[0]:
            doomed = [r.seqid.num for r in ch
                      if not (lo <= r.seqid.num + FIRST_RESIDUE - 1 <= hi)]
            for num in reversed(doomed):
                for i, r in enumerate(ch):
                    if r.seqid.num == num:
                        del ch[i]
                        break

    st.setup_entities()
    st.write_pdb(str(dest))


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--min-iptm", type=float, default=0.60)
    ap.add_argument("--keep", nargs=2, type=int, metavar=("LO", "HI"),
                    default=None,
                    help="restrict every chain to this residue range before "
                         "computing accessibility, e.g. --keep 44 110 to remove "
                         "the flexible N-terminal region. Use as a control: if "
                         "a residue's accessibility rises when the tails are "
                         "removed, they were occluding it.")
    args = ap.parse_args()
    keep = tuple(args.keep) if args.keep else None

    root = Path(args.results_dir).expanduser().resolve()
    out = (Path(args.out).expanduser().resolve() if args.out
           else root / "sasa_ensemble")
    tag = f"_keep{args.keep[0]}-{args.keep[1]}" if args.keep else ""
    tmpdir = out / f"_relabelled_pdb{tag}"
    tmpdir.mkdir(parents=True, exist_ok=True)

    backend, backend_name = pick_backend()
    print(f"Reading : {root}")
    print(f"Writing : {out}")
    print(f"SASA    : {backend_name}")
    if keep:
        print(f"Chains  : truncated to residues {keep[0]}-{keep[1]} "
              f"(occlusion control)")
    else:
        print("Chains  : full length")
    print()

    wanted = sorted({r for _, rs, _ in GROUPS for r in rs})

    rows = []          # one per model per chain per residue
    n_total = n_ap = n_accepted = 0
    excluded = []

    for cif in sorted(root.glob("results_*_seed*/**/*model_*.cif")):
        n_total += 1
        seed = next((int(p.split("_seed")[-1]) for p in cif.parts
                     if "_seed" in p), None)
        sample = cif.stem.split("_model_")[-1]

        conf_path = cif.parent / ("confidence_" + cif.name.replace(".cif", ".json"))
        iptm = (json.loads(conf_path.read_text()).get("iptm")
                if conf_path.exists() else None)

        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        cls, mapping = classify_and_map(st)
        if cls == "AP":
            n_ap += 1

        if cls != "AP":
            excluded.append((seed, sample, f"topology {cls}"))
            continue
        if iptm is None or iptm < args.min_iptm:
            excluded.append((seed, sample, f"ipTM {iptm:.3f}"))
            continue
        n_accepted += 1

        pdb = tmpdir / f"seed{seed:02d}_m{sample}.pdb"
        relabel_and_write(st, mapping, pdb, keep)

        areas = backend(pdb)
        st2 = gemmi.read_structure(str(pdb))
        st2.setup_entities()
        for ch in st2[0]:
            for res in ch:
                num = res.seqid.num + FIRST_RESIDUE - 1
                if num not in wanted:
                    continue
                a = areas.get((ch.name, res.seqid.num))
                if a is None:
                    continue
                mx = MAX_ASA.get(res.name)
                rows.append({"seed": seed, "sample": sample, "chain": ch.name,
                             "residue": num, "resn": res.name,
                             "abs_sasa": round(a, 2),
                             "rel_sasa": round(100 * a / mx, 1) if mx else None})

    if not rows:
        sys.exit("No accepted models found.")

    print(f"{n_total} models, {n_ap} antiparallel, {n_accepted} accepted "
          f"(AP and ipTM >= {args.min_iptm})")
    if excluded:
        print("  excluded:")
        for seed, sample, why in excluded:
            print(f"    seed {seed} sample {sample}: {why}")

    # ---- summarise
    def summarise(sel):
        v = np.array([r["rel_sasa"] for r in sel if r["rel_sasa"] is not None])
        return v.mean(), v.std(), v.min(), v.max()

    summary = []
    for title, resnums, note in GROUPS:
        print(f"\n{title}")
        print(f"  ({note})")
        print(f"    {'residue':>10s}{'chain':>7s}{'mean':>8s}{'SD':>7s}"
              f"{'min':>7s}{'max':>7s}   verdict")
        for num in resnums:
            for ch in ("A", "B", "C", "all"):
                sel = [r for r in rows if r["residue"] == num
                       and (ch == "all" or r["chain"] == ch)]
                if not sel:
                    continue
                m, sd, lo, hi = summarise(sel)
                resn = sel[0]["resn"]
                verdict = ("buried" if m < BURIED_BELOW else
                           "exposed" if m > EXPOSED_ABOVE else "intermediate")
                flag = "  <-- spread" if (hi - lo) > 30 else ""
                label = f"{resn}{num}" if ch == "A" else ""
                print(f"    {label:>10s}{ch:>7s}{m:7.1f}%{sd:6.1f}"
                      f"{lo:7.1f}{hi:7.1f}   {verdict}{flag}")
                summary.append({"group": title, "residue": num, "resn": resn,
                                "chain": ch, "mean_rel": round(m, 1),
                                "sd": round(sd, 1), "min": round(lo, 1),
                                "max": round(hi, 1), "verdict": verdict,
                                "n": len(sel)})
            print()

    # ---- write
    with open(out / f"sasa_per_model{tag}.tsv", "w") as fh:
        keys = ["seed", "sample", "chain", "residue", "resn",
                "abs_sasa", "rel_sasa"]
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[k]) for k in keys) + "\n")

    with open(out / f"sasa_summary{tag}.tsv", "w") as fh:
        keys = ["group", "residue", "resn", "chain", "mean_rel", "sd",
                "min", "max", "verdict", "n"]
        fh.write("\t".join(keys) + "\n")
        for r in summary:
            fh.write("\t".join(str(r[k]) for k in keys) + "\n")

    try:
        import pandas as pd
        meta = pd.DataFrame([
            ("Models found", n_total),
            ("Antiparallel", n_ap),
            ("Accepted", n_accepted),
            ("Acceptance", f"antiparallel and ipTM >= {args.min_iptm}"),
            ("SASA method", backend_name),
            ("Chain extent",
             f"residues {keep[0]}-{keep[1]}" if keep else "full length"),
            ("Relative to", "Tien et al. (2013) theoretical maxima"),
            ("Buried threshold", f"< {BURIED_BELOW}% relative"),
            ("Exposed threshold", f"> {EXPOSED_ABOVE}% relative"),
        ], columns=["Item", "Value"])
        with pd.ExcelWriter(out / f"sasa{tag}.xlsx", engine="openpyxl") as xw:
            meta.to_excel(xw, sheet_name="Method", index=False)
            pd.DataFrame(summary).to_excel(xw, sheet_name="Summary", index=False)
            pd.DataFrame(rows).to_excel(xw, sheet_name="Per model", index=False)
            for ws in xw.book.worksheets:
                ws.freeze_panes = "A2"
        print(f"\n  Excel: {out / f'sasa{tag}.xlsx'}")
    except ImportError:
        pass

    print(f"\nWritten to {out}/")
    print(f"  sasa_per_model{tag}.tsv, sasa_summary{tag}.tsv, "
          f"sasa{tag}.xlsx, _relabelled_pdb{tag}/")
    print("\nA large min-max spread for a residue means its accessibility is "
          "\n  model-dependent - expected where the ensemble RMSF is high, and "
          "\n  a reason not to quote a single model's value for that position.")


if __name__ == "__main__":
    main()
