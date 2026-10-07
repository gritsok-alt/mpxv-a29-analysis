#!/usr/bin/env python3
"""
Substitution-saturation test for a protein alignment and its ML tree.

Plots observed pairwise p-distance against pairwise patristic (tree) distance
for every taxon pair, and quantifies the saturation curve:

  * initial slope      least-squares fit over pairs below --linear-max
  * plateau slope      least-squares fit over pairs at or above --plateau-min
  * plateau level      mean p-distance of pairs at or above --plateau-min
  * plateau occupancy  fraction of pairs at or above --plateau-min
  * random expectation 1 - sum(f_i^2) over the amino-acid composition of the
                       alignment, i.e. the p-distance expected between two
                       sequences drawn independently from that composition

Reproduces Fig. S8 of the MPXV A29 (OPG154) orthologue analysis.

Usage
-----
  python saturation_analysis.py --msa A29L_msa_linsi.fasta --tree A29L.treefile

Outputs (written to --outdir, default ./saturation_out)
  figS8_saturation.png          the figure
  saturation_pairwise.csv       one row per taxon pair (3 cols)
  saturation_stats.csv          every reported statistic, one per row
  saturation_stats.json         same, machine-readable
"""

import argparse
import itertools
import json
import os
from collections import Counter

import numpy as np

# ---------------------------------------------------------------- I/O


def read_fasta(path):
    """Return {header: sequence} preserving file order. Headers are the full
    line after '>' with surrounding whitespace stripped."""
    seqs, h = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                h = line[1:].strip()
                if h in seqs:
                    raise ValueError(f"duplicate FASTA header: {h!r}")
                seqs[h] = []
            elif h is None:
                raise ValueError("sequence data before first '>' header")
            else:
                seqs[h].append(line)
    out = {k: "".join(v).upper() for k, v in seqs.items()}
    lens = {len(s) for s in out.values()}
    if len(lens) != 1:
        raise ValueError(
            f"sequences are not all the same length ({sorted(lens)}); "
            "input must be an ALIGNED FASTA"
        )
    return out


class Node:
    __slots__ = ("ch", "name", "bl", "parent")

    def __init__(self):
        self.ch = []
        self.name = None
        self.bl = 0.0
        self.parent = None


def read_newick(path):
    """Parse a Newick tree with branch lengths. Internal node labels (support
    values, as written by IQ-TREE) are read and discarded. Quoted tip labels
    are not supported; use safe tip names (see README)."""
    s = open(path).read().strip()
    if not s.endswith(";"):
        s += ";"
    i = 0

    def node():
        nonlocal i
        n = Node()
        if s[i] == "(":
            i += 1
            while True:
                n.ch.append(node())
                if s[i] == ",":
                    i += 1
                else:
                    break
            if s[i] != ")":
                raise ValueError(f"expected ')' at offset {i}")
            i += 1
        j = i
        while i < len(s) and s[i] not in "(),:;":
            i += 1
        label = s[j:i]
        if not n.ch:
            n.name = label
        if i < len(s) and s[i] == ":":
            i += 1
            j = i
            while i < len(s) and s[i] not in "(),;":
                i += 1
            n.bl = float(s[j:i])
        return n

    root = node()

    def link(n, p=None):
        n.parent = p
        for c in n.ch:
            link(c, n)

    link(root)
    return root


def tips(root):
    out = {}

    def walk(n):
        if not n.ch:
            if n.name in out:
                raise ValueError(f"duplicate tip label: {n.name!r}")
            out[n.name] = n
        for c in n.ch:
            walk(c)

    walk(root)
    return out


# ---------------------------------------------------------------- distances


def patristic_matrix(tip_nodes, names):
    """Pairwise patristic distance = sum of branch lengths on the path between
    two tips, computed as d(a,mrca) + d(b,mrca) where the MRCA is the shallowest
    node on both root-paths. Root placement does not affect the result."""
    paths = {}
    for nm in names:
        n, d, acc = tip_nodes[nm], 0.0, {}
        while n is not None:
            acc[id(n)] = d
            d += n.bl
            n = n.parent
        paths[nm] = acc
    out = {}
    for a, b in itertools.combinations(names, 2):
        Pa, Pb = paths[a], paths[b]
        out[(a, b)] = min(Pa[k] + Pb[k] for k in Pa if k in Pb)
    return out


def p_distance_matrix(seqs, names, gap_chars="-.~"):
    """Pairwise p-distance with PAIRWISE gap deletion: for each pair, compare
    only the columns where BOTH sequences carry a residue, and return the
    fraction of those columns that differ. Columns gapped in either sequence
    are ignored, so the denominator varies between pairs."""
    arrs = {nm: np.frombuffer(seqs[nm].encode(), dtype="S1") for nm in names}
    gaps = {g.encode() for g in gap_chars}
    present = {nm: ~np.isin(arrs[nm], list(gaps)) for nm in names}
    out, ncomp = {}, {}
    for a, b in itertools.combinations(names, 2):
        m = present[a] & present[b]
        n = int(m.sum())
        if n == 0:
            raise ValueError(f"no shared ungapped columns for {a!r} and {b!r}")
        out[(a, b)] = float((arrs[a][m] != arrs[b][m]).sum()) / n
        ncomp[(a, b)] = n
    return out, ncomp


def random_expectation(seqs, gap_chars="-.~"):
    """p-distance expected between two sequences drawn independently from the
    observed amino-acid composition of the alignment: 1 - sum_i f_i^2."""
    freq = Counter()
    for s in seqs.values():
        freq.update(c for c in s if c not in gap_chars)
    tot = sum(freq.values())
    f = np.array([v / tot for v in freq.values()])
    return float(1.0 - (f ** 2).sum())


# ---------------------------------------------------------------- stats


def lsq_slope(x, y):
    """Ordinary least-squares slope and intercept."""
    if len(x) < 3:
        return float("nan"), float("nan")
    b, a = np.polyfit(x, y, 1)
    return float(b), float(a)


def compute_stats(x, y, rand_p, linear_max, plateau_min):
    s_lin, i_lin = lsq_slope(x[x < linear_max], y[x < linear_max])
    s_plat, i_plat = lsq_slope(x[x >= plateau_min], y[x >= plateau_min])
    st = {
        "n_pairs": int(len(x)),
        "random_expectation_p": rand_p,
        "linear_window_max": linear_max,
        "plateau_window_min": plateau_min,
        "n_pairs_linear_window": int((x < linear_max).sum()),
        "initial_slope": s_lin,
        "initial_intercept": i_lin,
        "n_pairs_plateau": int((x >= plateau_min).sum()),
        "plateau_fraction_of_pairs": float((x >= plateau_min).mean()),
        "plateau_mean_p": float(y[x >= plateau_min].mean()),
        "plateau_slope": s_plat,
        "plateau_intercept": i_plat,
        "plateau_mean_p_over_random": float(y[x >= plateau_min].mean() / rand_p),
        "max_patristic": float(x.max()),
        "max_p_distance": float(y.max()),
        "mean_p_all_pairs": float(y.mean()),
    }
    # Robustness: the same quantities at alternative window edges. These are
    # the numbers a reader may compute with a different threshold in mind; the
    # reported value must match the threshold stated in the text.
    for cut in (1.0, 2.0, 4.0):
        m = x >= cut
        if m.sum() >= 3:
            st[f"alt_fraction_pairs_ge_{cut:g}"] = float(m.mean())
            st[f"alt_mean_p_ge_{cut:g}"] = float(y[m].mean())
        m2 = x < cut
        if m2.sum() >= 3:
            st[f"alt_slope_below_{cut:g}"] = lsq_slope(x[m2], y[m2])[0]
            st[f"alt_mean_p_below_{cut:g}"] = float(y[m2].mean())
    return st


def binned_means(x, y, bins, min_per_bin=10):
    ctr, mn, nb = [], [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (x >= lo) & (x < hi)
        if m.sum() >= min_per_bin:
            ctr.append(float(x[m].mean()))
            mn.append(float(y[m].mean()))
            nb.append(int(m.sum()))
    return np.array(ctr), np.array(mn), nb


# ---------------------------------------------------------------- figure


import numpy as np

GREY = "#888888"


def apply_figure_style(sizes=(9, 8, 8)):
    import matplotlib as mpl

    base, secondary, tick = sizes
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Tinos", "Liberation Serif",
                           "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": base,
            "axes.labelsize": base,
            "legend.fontsize": secondary,
            "xtick.labelsize": tick,
            "ytick.labelsize": tick,
            "axes.linewidth": 0.7,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "legend.frameon": False,
            "figure.dpi": 200,
            "savefig.dpi": 1000,
            "savefig.bbox": "tight",
            "lines.linewidth": 1.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def _eq(slope, intercept, dp=2):
    sign = "+" if intercept >= 0 else "\u2212"
    return f"$y = {slope:.{dp}f}x\\ {sign}\\ {abs(intercept):.2f}$"


def make_figure(x, y, st, bins, outpath, title=None, dpi=1000):
    import matplotlib.pyplot as plt
    from pathlib import Path

    apply_figure_style()
    ctr, mn, _ = binned_means(x, y, bins)  # noqa: F821  (from the module)
    xmax = float(np.ceil(x.max() * 1.03))
    rand_p = st["random_expectation_p"]
    pmin = st["plateau_window_min"]
    lmax = st["linear_window_max"]

    fig, ax = plt.subplots(figsize=(4.8, 3.6))

    # Pairwise comparisons as a density cloud: hexagonal binning, so the
    # concentration of pairs is visible instead of 1,770 overplotted markers.
    hb = ax.hexbin(x, y, gridsize=46, cmap="Greys", mincnt=1,
                   linewidths=0, zorder=1)
    hb.set_alpha(0.85)

    # Compositional ceiling.
    ax.axhline(rand_p, ls=(0, (6, 3)), lw=1.0, color="black", zorder=4)
    ax.text(xmax, rand_p + 0.018,
            f"compositional ceiling ({rand_p:.2f})",
            ha="right", va="bottom", fontsize=7.5, color="black")

    # Boundary of the plateau window.
    ax.axvline(pmin, ls=(0, (1, 3)), lw=0.8, color=GREY, zorder=2)

    # Fitted relationships, black, distinguished by dash pattern.
    xx = np.linspace(0, lmax * 1.75, 50)
    ax.plot(xx, st["initial_slope"] * xx + st["initial_intercept"],
            ls=(0, (1, 2)), lw=1.4, color="black", zorder=5)
    xp = np.linspace(pmin, xmax, 50)
    ax.plot(xp, st["plateau_slope"] * xp + st["plateau_intercept"],
            ls=(0, (5, 2)), lw=1.4, color="black", zorder=5)

    # Binned means.
    ax.plot(ctr, mn, "-", color="black", lw=1.0, zorder=6)
    ax.plot(ctr, mn, "o", color="black", ms=3.4, zorder=6)

    # Equations, placed in open space.
    ax.text(0.30, 0.86, _eq(st["initial_slope"], st["initial_intercept"]),
            transform=ax.transAxes, fontsize=8, ha="left", va="center")
    ax.text(0.30, 0.80,
            f"below {lmax:g} substitutions per site",
            transform=ax.transAxes, fontsize=7, ha="left", va="center",
            color="#333333")

    ax.text(0.30, 0.70, _eq(st["plateau_slope"], st["plateau_intercept"], dp=3),
            transform=ax.transAxes, fontsize=8, ha="left", va="center")
    ax.text(0.30, 0.64,
            f"at or above {pmin:g} substitutions per site",
            transform=ax.transAxes, fontsize=7, ha="left", va="center",
            color="#333333")

    ax.set_xlabel("Patristic distance (substitutions per site)")
    ax.set_ylabel("Observed p-distance")
    ax.set_ylim(0, 1.0)
    ax.set_xlim(-0.03 * xmax, xmax * 1.03)
    fig.tight_layout()

    out = Path(outpath)
    fig.savefig(out, dpi=dpi)
    for ext in (".pdf", ".svg"):
        fig.savefig(out.with_suffix(ext))

    r = fig.canvas.get_renderer()
    boxes = [(t, t.get_window_extent(r)) for a in fig.axes for t in a.texts
             if t.get_text().strip()]
    clashes = [(a.get_text()[:18], b.get_text()[:18])
               for i, (a, ba) in enumerate(boxes) for b, bb in boxes[i + 1:]
               if ba.overlaps(bb)]
    return fig, clashes


# ---------------------------------------------------------------- main


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--msa", required=True, help="aligned FASTA (all sequences equal length)")
    ap.add_argument("--tree", required=True, help="Newick tree with branch lengths")
    ap.add_argument("--outdir", default="saturation_out")
    ap.add_argument("--linear-max", type=float, default=1.0,
                    help="upper edge of the initial-slope fitting window "
                         "(substitutions per site; default 1.0)")
    ap.add_argument("--plateau-min", type=float, default=2.0,
                    help="lower edge of the plateau window (default 2.0)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--figname", default="figS8_saturation.png")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    seqs = read_fasta(args.msa)
    tip_nodes = tips(read_newick(args.tree))

    # Tip labels must match FASTA headers exactly. IQ-TREE rewrites characters
    # it considers unsafe, so this is the usual failure point.
    only_msa = sorted(set(seqs) - set(tip_nodes))
    only_tree = sorted(set(tip_nodes) - set(seqs))
    if only_msa or only_tree:
        raise SystemExit(
            "tip labels do not match FASTA headers.\n"
            f"  {len(only_msa)} only in MSA:  {only_msa[:4]}\n"
            f"  {len(only_tree)} only in tree: {only_tree[:4]}\n"
            "Rename both to a common set of safe identifiers (see README)."
        )

    names = list(seqs)
    pat = patristic_matrix(tip_nodes, names)
    pdist, ncomp = p_distance_matrix(seqs, names)

    pairs = list(itertools.combinations(names, 2))
    x = np.array([pat[p] for p in pairs])
    y = np.array([pdist[p] for p in pairs])
    rand_p = random_expectation(seqs)

    st = compute_stats(x, y, rand_p, args.linear_max, args.plateau_min)
    st["n_sequences"] = len(names)
    st["n_alignment_columns"] = len(next(iter(seqs.values())))
    st["msa_file"] = os.path.basename(args.msa)
    st["tree_file"] = os.path.basename(args.tree)

    bins = np.array([0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10])
    if x.max() > bins[-1]:
        bins = np.append(bins, np.ceil(x.max()))

    pair_csv = os.path.join(args.outdir, "saturation_pairwise.csv")
    with open(pair_csv, "w") as fh:
        fh.write("taxon_a,taxon_b,patristic,p_distance,n_compared_columns\n")
        for p in pairs:
            fh.write(f'"{p[0]}","{p[1]}",{pat[p]:.6f},{pdist[p]:.6f},{ncomp[p]}\n')

    with open(os.path.join(args.outdir, "saturation_stats.json"), "w") as fh:
        json.dump(st, fh, indent=1, sort_keys=True)
    with open(os.path.join(args.outdir, "saturation_stats.csv"), "w") as fh:
        fh.write("statistic,value\n")
        for k in sorted(st):
            fh.write(f"{k},{st[k]}\n")

    figpath = os.path.join(args.outdir, args.figname)
    _, clashes = make_figure(x, y, st, bins, figpath, title=args.title)

    print(f"{st['n_sequences']} sequences, {st['n_alignment_columns']} columns, "
          f"{st['n_pairs']} pairs")
    for k in ("random_expectation_p", "initial_slope", "initial_intercept",
              "plateau_mean_p", "plateau_slope", "plateau_fraction_of_pairs",
              "plateau_mean_p_over_random", "max_patristic"):
        print(f"  {k:28s} {st[k]:.4f}")
    if clashes:
        print(f"WARNING: overlapping figure text: {clashes}")
    print(f"wrote {figpath}, saturation_pairwise.csv, saturation_stats.{{csv,json}}")
    return st


if __name__ == "__main__":
    main()
