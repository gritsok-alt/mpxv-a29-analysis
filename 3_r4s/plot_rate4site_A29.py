#!/usr/bin/env python3
"""
Reproduce rate4site_A29_MPXV.png (v2) — site-specific evolutionary rates of the
MPXV A29 protein, referenced to YFT67316.1 (110 aa).

DATA SOURCE
-----------
The figure is drawn entirely from the primary rate4site run:

    r4s_run/r60_jtt.res      (inside rate4site_A29_run.tar.gz)

which was produced by:

    rate4site -s aln.fasta -t tree.nwk -a T60 -Mj -ib -k 16 \
              -o r60_jtt.res -y r60_jtt.unnorm -x r60_jtt.tree

    -a T60  reference = Orthopoxvirus|YFT67316.1|Orthopoxvirus_monkeypox
    -Mj     JTT replacement matrix
    -ib     empirical-Bayes rate inference
    -k 16   16 discrete gamma categories

Every plotted quantity (score, 25-75% posterior interval, grade, low-confidence
flag) is derived from that one file. The WAG / fixed-branch-length / Q.INSECT
columns in the published CSV are sensitivity checks and are NOT used here.

Two equivalent inputs are accepted:
    --res  r60_jtt.res                      (raw rate4site output; grades re-derived)
    --csv  rate4site_A29_MPXV_scores.csv    (published table; grades read directly)
Both give a byte-identical figure.

USAGE
-----
    python plot_rate4site_A29.py --res r60_jtt.res
    python plot_rate4site_A29.py --csv rate4site_A29_MPXV_scores.csv
    python plot_rate4site_A29.py --res r60_jtt.res --no-font-fetch   # use local Times

Outputs rate4site_A29_MPXV.png (1000 dpi) and rate4site_A29_MPXV.pdf.
"""

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D

# --------------------------------------------------------------------------
# 1. Times New Roman
# --------------------------------------------------------------------------
# On a machine with real Times New Roman installed this is a no-op. In a bare
# Linux container nothing Times-metric exists, so we fall back to Tinos
# (Apache-2.0, metrically compatible with Times New Roman) and rewrite its
# internal name records so matplotlib's family lookup resolves.

TINOS_URL = "https://raw.githubusercontent.com/googlefonts/tinos/main/fonts/ttf/Tinos-{}.ttf"
FACES = {"Regular": "", "Bold": "Bold", "Italic": "Italic", "BoldItalic": "BoldItalic"}


def times_available():
    try:
        p = fm.findfont("Times New Roman", fallback_to_default=False)
        return bool(p) and "DejaVu" not in Path(p).name
    except Exception:
        return False


def setup_times_font(fetch=True, fontdir="fonts/tnr"):
    """Ensure the family name 'Times New Roman' resolves. Returns True on success."""
    if times_available():
        return True
    fd = Path(fontdir)
    if fd.is_dir() and list(fd.glob("*.ttf")):
        for f in fd.glob("*.ttf"):
            fm.fontManager.addfont(str(f))
        return times_available()
    if not fetch:
        return False

    import urllib.request
    from fontTools.ttLib import TTFont

    fd.mkdir(parents=True, exist_ok=True)
    raw = Path("fonts/tinos")
    raw.mkdir(parents=True, exist_ok=True)
    for face, suffix in FACES.items():
        src = raw / f"Tinos-{face}.ttf"
        if not src.exists():
            urllib.request.urlretrieve(TINOS_URL.format(face), src)
        # Rewrite name records: family (1), full name (4), PostScript (6),
        # typographic family (16) if present.
        ft = TTFont(str(src))
        family = "Times New Roman"
        full = (family + " " + suffix).strip()
        ps = "TimesNewRoman" + suffix
        for rec in ft["name"].names:
            if rec.nameID == 1:
                rec.string = family
            elif rec.nameID == 4:
                rec.string = full
            elif rec.nameID == 6:
                rec.string = ps
            elif rec.nameID == 16:
                rec.string = family
        out = fd / f"TimesNewRoman{suffix}.ttf"
        ft.save(str(out))
        ft.close()
        fm.fontManager.addfont(str(out))
    fm._load_fontmanager(try_read_cache=False)
    for f in fd.glob("*.ttf"):
        fm.fontManager.addfont(str(f))
    return times_available()


# --------------------------------------------------------------------------
# 2. Figure style  (from the figure-style skill, inlined so this stands alone)
# --------------------------------------------------------------------------
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    base, secondary, tick = sizes
    boxed = frame == "boxed"
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, fontsize=None):
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    ax.text(dx, dy, letter.lower(), transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


# --------------------------------------------------------------------------
# 3. Data
# --------------------------------------------------------------------------
# rate4site .res columns: pos  aa  score  [qq_lo, qq_hi]  std  n_nongap/n_total
R4S_RE = re.compile(
    r"\s*(\d+)\s+(\S)\s+(-?[\d.eE+-]+)\s+"
    r"\[\s*(-?[\d.eE+-]+),\s*(-?[\d.eE+-]+)\]\s+"
    r"(-?[\d.eE+-]+)\s+(\d+)/(\d+)"
)


def read_r4s(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        m = R4S_RE.match(line)
        if m:
            rows.append(dict(
                pos=int(m.group(1)), aa=m.group(2), score=float(m.group(3)),
                qq_lo=float(m.group(4)), qq_hi=float(m.group(5)),
                std=float(m.group(6)),
                msa_n=int(m.group(7)), msa_tot=int(m.group(8)),
            ))
    return pd.DataFrame(rows)


def add_grades(d):
    """ConSurf-style 9-bin discretization of the normalized rate.

    rate4site's score is normalized to mean 0, sd 1 across the reference sites,
    and LOWER = MORE CONSERVED. Bins are equal-interval over the observed score
    range; grade 9 = most conserved. A site is flagged low-confidence when its
    25-75% posterior interval spans more than a third of the score range.
    """
    edges = np.linspace(d.score.min(), d.score.max(), 10)
    d["grade"] = 9 - np.clip(np.digitize(d.score, edges[1:-1]), 0, 8)
    span = d.score.max() - d.score.min()
    d["unreliable"] = (d.qq_hi - d.qq_lo) > span / 3
    return d


def load(args):
    if args.res:
        d = add_grades(read_r4s(args.res))
    else:
        c = pd.read_csv(args.csv)
        d = pd.DataFrame(dict(
            pos=c.pos_ref_MPXV_A29, aa=c.aa, score=c.rate_score,
            qq_lo=c.qq25, qq_hi=c.qq75,
            grade=c.consurf_grade.astype(int),
            unreliable=c.low_confidence.astype(bool),
        ))
    assert len(d) == 110, f"expected 110 reference sites, got {len(d)}"
    assert d.pos.tolist() == list(range(1, 111)), "positions must be 1..110"
    return d


# --------------------------------------------------------------------------
# 4. Figure
# --------------------------------------------------------------------------
# ConSurf palette: 1 = variable (cyan) ... 9 = conserved (maroon)
CONSURF_HEX = ["#10C8D1", "#8CFFFF", "#D7FFFF", "#EAFFFF", "#FFFFFF",
               "#FCEDF4", "#FAC0DB", "#F07DAB", "#A02560"]
GRADE_COLOR = {g: CONSURF_HEX[g - 1] for g in range(1, 10)}
UNRELIABLE = "#FFFF96"          # ConSurf "insufficient data" yellow
PERROW = 37                      # 110 residues over 3 lines


def make_figure(d, have_times):
    font = "Times New Roman" if have_times else None
    apply_figure_style(font=font, sizes=(11, 9.5, 8.5))
    if have_times:
        mpl.rcParams["font.serif"] = ["Times New Roman", "Tinos", "DejaVu Serif"]
    mpl.rcParams["mathtext.fontset"] = "stix"

    nrows = int(np.ceil(len(d) / PERROW))
    bar_col = [UNRELIABLE if u else GRADE_COLOR[g]
               for g, u in zip(d.grade, d.unreliable)]

    fig = plt.figure(figsize=(7.4, 5.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.72], hspace=0.42)

    # ---- Panel a: one column per residue, conserved pointing down ----------
    axA = fig.add_subplot(gs[0])
    x = d.pos.values
    axA.vlines(x, d.qq_lo, d.qq_hi, color="#9AA5B1", lw=0.85, zorder=1)
    axA.bar(x, d.score, width=0.78, color=bar_col,
            edgecolor="#3C4650", lw=0.35, zorder=2)
    axA.axhline(0, color="#3C4650", lw=0.8, zorder=3)
    axA.set_xlim(0.0, 111.0)
    axA.set_ylim(d.qq_lo.min() - 0.18, d.qq_hi.max() + 0.18)
    axA.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110])
    axA.set_xlabel("Residue position in MPXV A29 (YFT67316.1, 110 aa)")
    axA.set_ylabel("Normalized evolutionary rate")
    axA.set_title("Site-specific rates across A29: a variable N-terminal half "
                  "and a conserved C-terminal half", loc="left", pad=22)
    axA.text(0.155, 0.965, "more variable  \u2191", transform=axA.transAxes,
             ha="left", va="top", fontsize=8.5, color=META_GREY)
    axA.text(0.155, 0.035, "more conserved  \u2193", transform=axA.transAxes,
             ha="left", va="bottom", fontsize=8.5, color=META_GREY)
    axA.legend(handles=[
        mpatches.Patch(facecolor="#F3F5F7", edgecolor="#3C4650", lw=0.35,
                       label="Rate per residue (fill = grade, panel b)"),
        Line2D([], [], color="#9AA5B1", lw=0.85,
               label="25\u201375 % posterior interval")],
        loc="lower left", bbox_to_anchor=(0.0, 1.002), ncol=2, frameon=False,
        handlelength=1.3, columnspacing=1.4, borderaxespad=0.0,
        handletextpad=0.5)
    panel_letter(axA, "a")

    # ---- Panel b: graded sequence block, 3 lines of 37 -------------------
    axB = fig.add_subplot(gs[1])
    axB.set_axis_off()
    cw, rh, rgap = 1.0, 1.0, 1.85
    for r in d.itertuples():
        ri, ci = divmod(r.pos - 1, PERROW)
        col = UNRELIABLE if r.unreliable else GRADE_COLOR[r.grade]
        y0 = -ri * rgap
        axB.add_patch(mpatches.Rectangle((ci * cw, y0), cw, rh,
                                         facecolor=col, edgecolor="#5A6570", lw=0.3))
        axB.text(ci * cw + cw / 2, y0 + rh / 2, r.aa, ha="center", va="center",
                 fontsize=5.6,
                 color="white" if (r.grade == 9 and not r.unreliable) else "black")
        axB.text(ci * cw + cw / 2, y0 + rh + 0.10, str(r.grade),
                 ha="center", va="bottom", fontsize=4.4, color="#3C4650")
    for ri in range(nrows):
        s0 = ri * PERROW + 1
        e0 = min(s0 + PERROW - 1, len(d))
        y0 = -ri * rgap
        axB.text(-0.3, y0 + rh / 2, str(s0), ha="right", va="center", fontsize=6.0)
        axB.text(PERROW * cw + 0.3, y0 + rh / 2, str(e0), ha="left", va="center",
                 fontsize=6.0)
    axB.set_xlim(-2.4, PERROW * cw + 2.4)
    axB.set_ylim(-(nrows - 1) * rgap - 1.75, rh + 1.15)
    axB.set_title("Conservation grade per residue "
                  "(9 = most conserved, 1 = most variable)",
                  loc="left", x=-0.004, y=1.0)

    # grade key + low-confidence swatch
    ky = -(nrows - 1) * rgap - 1.05
    kw = 0.95
    for g in range(1, 10):
        axB.add_patch(mpatches.Rectangle(((g - 1) * kw, ky), kw, 0.5,
                                         facecolor=GRADE_COLOR[g],
                                         edgecolor="#5A6570", lw=0.3))
        axB.text((g - 1) * kw + kw / 2, ky + 0.25, str(g), ha="center",
                 va="center", fontsize=5.6, color="white" if g == 9 else "black")
    axB.text(-0.25, ky + 0.25, "variable", ha="right", va="center", fontsize=6.5)
    axB.text(9 * kw + 0.25, ky + 0.25, "conserved", ha="left", va="center",
             fontsize=6.5)
    axB.add_patch(mpatches.Rectangle((16.0, ky), kw, 0.5, facecolor=UNRELIABLE,
                                     edgecolor="#5A6570", lw=0.3))
    axB.text(16.0 + kw + 0.3, ky + 0.25,
             "insufficient data (wide posterior interval)",
             ha="left", va="center", fontsize=6.5)
    panel_letter(axB, "b")

    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--res", help="raw rate4site output, e.g. r60_jtt.res")
    src.add_argument("--csv", help="published table rate4site_A29_MPXV_scores.csv")
    ap.add_argument("--out", default="rate4site_A29_MPXV",
                    help="output basename (default: %(default)s)")
    ap.add_argument("--dpi", type=int, default=1000)
    ap.add_argument("--no-font-fetch", action="store_true",
                    help="do not download Tinos; use a locally installed Times")
    args = ap.parse_args()

    have_times = setup_times_font(fetch=not args.no_font_fetch)
    if not have_times:
        print("warning: 'Times New Roman' did not resolve; using the matplotlib "
              "default sans-serif. Install Times New Roman or drop "
              "--no-font-fetch.", file=sys.stderr)

    d = load(args)
    fig = make_figure(d, have_times)
    fig.savefig(f"{args.out}.png", dpi=args.dpi, bbox_inches="tight")
    fig.savefig(f"{args.out}.pdf", bbox_inches="tight")
    print(f"wrote {args.out}.png and {args.out}.pdf "
          f"({len(d)} sites; {int(d.unreliable.sum())} flagged low-confidence)")


if __name__ == "__main__":
    main()
