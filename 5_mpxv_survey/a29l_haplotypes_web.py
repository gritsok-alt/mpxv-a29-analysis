#!/usr/bin/env python3
"""
A29L haplotype survey from an NCBI Virus web export.

Expects headers of the form:
    >YDZ51540.1 |IMV surface fusion protein [Monkeypox virus]|IIb|2026-03-24
                                                              ^clade ^date

QC is applied as REPORTING then exclusion: every count is printed before
anything is removed, and every excluded sequence is written to disk.

Usage:  python a29l_haplotypes_web.py [sequences.fasta] [--outdir results]
"""
import sys, re, collections, argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import SeqIO

# --------------------------------------------------------------------------- #
REF = ("MDGTLFPGDDDLAIPATEFFSTKAAKNPETKREAIVKAYGDDNEETLKQRLTNLEKKITN"
       "ITTKFEQIEKCCKHNDEVLFRLENHAETLRAAMISLAKKIDVQTGRRPYE")      # Q77HM6, clade Ia
STD_AA   = set("ACDEFGHIKLMNPQRSTVWY")
VALID    = {"I", "IA", "IB", "II", "IIA", "IIB"}
MIN_ID   = 0.90
# A29 domain architecture (Wang et al. 2014; Chang et al. 2013)
DOMAINS  = [(1, 20, "NTR", "#E8E8E8"), (21, 34, "HBD", "#FDBE85"), (35, 42, "sp", "#E8E8E8"),
            (43, 84, "CCD", "#B7E2B1"), (85, 110, "LZD", "#C6DBEF")]
OKABE    = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#8C6D31"]

ap = argparse.ArgumentParser()
ap.add_argument("fasta", nargs="?", default="sequences.fasta")
ap.add_argument("--outdir", default="results_web")
A = ap.parse_args()
OUT = Path(A.outdir); OUT.mkdir(exist_ok=True)
log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

# ============================== 1 · PARSE ================================== #
rows = []
for rec in SeqIO.parse(A.fasta, "fasta"):
    p = rec.description.split("|")
    clade = p[2].strip() if len(p) > 2 else ""
    date  = p[3].strip() if len(p) > 3 else ""
    m = re.match(r"([A-Z]+)(\d+)", rec.id)                 # accession prefix + number
    rows.append({"accession": rec.id,
                 "prefix": m.group(1) if m else "",
                 "acc_num": int(m.group(2)) if m else np.nan,
                 "clade_raw": clade,
                 # 'Clade' occurs as a literal value - a field-label artefact of the export
                 "clade": clade if clade.upper() in VALID else "unassigned",
                 "date": date,
                 "year": int(date[:4]) if re.match(r"\d{4}", date) else np.nan,
                 "seq": str(rec.seq).upper()})
df = pd.DataFrame(rows)
log(f"[1] sequences read                       : {len(df):>7,}")

# ============================== 2 · QC ===================================== #
df["length"]   = df.seq.str.len()
df["bad_aa"]   = [sorted(set(s) - STD_AA) for s in df.seq]
df["identity"] = [sum(a == b for a, b in zip(s, REF)) / len(REF) if len(s) == len(REF) else 0.0
                  for s in df.seq]

log(f"[2] length distribution:")
for L, n in df.length.value_counts().sort_index().items():
    log(f"      {L:>4} aa : {n:>7,}{'   <- reference' if L == len(REF) else ''}")

n_len  = (df.length != len(REF)).sum()
n_amb  = (df.bad_aa.str.len() > 0).sum()
n_unrel= (df.identity < MIN_ID).sum()
log(f"[3] length != {len(REF)}                        : {n_len:>7,}")
log(f"[4] non-standard residues (X, B, Z, ...) : {n_amb:>7,}")
if n_amb:
    seen = collections.Counter(c for lst in df.bad_aa for c in lst)
    log(f"      characters seen: {dict(seen)}")
log(f"[5] identity < {MIN_ID:.0%} (unrelated protein)   : {n_unrel:>7,}")

df["exclude_reason"] = [
    ";".join(filter(None, [
        f"length_{r.length}" if r.length != len(REF) else "",
        f"nonstandard_{''.join(r.bad_aa)}" if r.bad_aa else "",
        "unrelated" if r.identity < MIN_ID else ""]))
    for r in df.itertuples()]
excluded = df[df.exclude_reason != ""]
d = df[df.exclude_reason == ""].reset_index(drop=True)
excluded.drop(columns=["seq"]).assign(sequence=excluded.seq).to_csv(
    OUT / "excluded_sequences.csv", index=False)
log(f"[6] FINAL curated set                    : {len(d):>7,}   "
    f"({len(excluded):,} excluded -> excluded_sequences.csv)")

log(f"\n[7] clade labels:")
for k, v in d.clade.value_counts().items():
    log(f"      {k:<12} {v:>7,}")

# ======================= 3 · VARIABLE POSITIONS ============================ #
M = np.array([list(s) for s in d.seq])
var = []
for i in range(len(REF)):
    c = collections.Counter(M[:, i])
    if len(c) > 1:
        maj, n = c.most_common(1)[0]
        var.append({"position": i + 1, "reference_aa": REF[i], "major_aa": maj, "major_n": n,
                    "minor_n": len(d) - n,
                    "minor": "; ".join(f"{a}={k}" for a, k in c.most_common()[1:]),
                    "domain": next((lab for s0, e0, lab, _ in DOMAINS if s0 <= i + 1 <= e0), "")})
var = pd.DataFrame(var)
var.to_csv(OUT / "variable_positions.csv", index=False)
log(f"\n[8] variable positions: {len(var)} of {len(REF)}")
log(var.drop(columns='minor_n').to_string(index=False))

# ============================ 4 · HAPLOTYPES =============================== #
counts = collections.Counter(d.seq)
dom = counts.most_common(1)[0][0]
hid = {s: f"H{i+1}" for i, (s, _) in enumerate(counts.most_common())}
d["haplotype"] = d.seq.map(hid)
d["res74"], d["res107"] = d.seq.str[73], d.seq.str[106]
d["pair"] = d.res74 + "74/" + d.res107 + "107"

haplo = pd.DataFrame([
    {"haplotype": hid[s], "n": n, "percent": round(100 * n / len(d), 3),
     "res74": s[73], "res107": s[106],
     "differences_vs_dominant": "; ".join(f"{dom[j]}{j+1}{s[j]}"
                                          for j in range(len(s)) if s[j] != dom[j]) or "DOMINANT"}
    for s, n in counts.most_common()])
haplo.to_csv(OUT / "haplotypes.csv", index=False)
log(f"\n[9] distinct haplotypes: {len(counts)}")
log(haplo.head(20).to_string(index=False))

# ========================= 5 · CLADE ASSOCIATION =========================== #
d["clade_group"] = d.clade.map(lambda c: "I" if c in ("I", "Ia", "Ib")
                               else ("II" if c in ("II", "IIa", "IIb") else "unassigned"))
ct_full = pd.crosstab(d.haplotype, d.clade).reindex(haplo.haplotype).fillna(0).astype(int)
ct_full.to_csv(OUT / "haplotype_by_clade.csv")

lab = d[d.clade_group != "unassigned"]
ct_pair = pd.crosstab(lab.clade_group, lab.pair)
log(f"\n[10] CLADE x residue-74/107 pair  (n = {len(lab):,} clade-assigned)")
log(ct_pair.to_string())
off = ct_pair.values.sum() - np.trace(ct_pair.values[np.ix_(
    range(len(ct_pair)), [list(ct_pair.columns).index(c) for c in ct_pair.idxmax(axis=1)])])
disc = ct_pair.values.sum() - sum(ct_pair.loc[r].max() for r in ct_pair.index)
log(f"     discordant sequences: {disc}")

un = d[d.clade_group == "unassigned"]
log(f"\n[11] unassigned ({len(un):,}) partition by pair: {un.pair.value_counts().to_dict()}")

# Because the pair->clade mapping is perfectly concordant in labelled data, clade can be
# imputed for unlabelled sequences. Used ONLY for denominators; always reported as imputed.
pair2clade = {p: ct_pair[p].idxmax() for p in ct_pair.columns}
d["clade_imputed"] = np.where(d.clade_group != "unassigned", d.clade_group,
                              d.pair.map(pair2clade).fillna("unassigned"))
log(f"     imputed via pair->clade {pair2clade}; "
    f"group sizes now {d.clade_imputed.value_counts().to_dict()}")

# Haplotypes that DEFINE a clade are not candidate 'emerging variants' - exclude them
# from the batch and temporal analyses, which exist to test low-frequency change.
clade_defining = {d[d.clade_imputed == g].haplotype.value_counts().idxmax()
                  for g in ("I", "II") if (d.clade_imputed == g).any()}
log(f"     clade-defining haplotypes excluded from [12]-[13]: {sorted(clade_defining)}")

# ================== 6 · BATCH-ARTEFACT DIAGNOSTIC ========================== #
# Submitter is not in the FASTA header, but consecutive accession numbers within
# a single prefix are a reliable proxy for one submission batch.
log("\n[12] SUBMISSION-BATCH DIAGNOSTIC for low-frequency haplotypes")
log("     A variant confined to one accession prefix and a narrow numeric span")
log("     is a single deposit, not evidence of spread.")
batch = []
cand = haplo.haplotype[(haplo.n >= 5) & (~haplo.haplotype.isin(clade_defining))]
for h in cand:
    sub = d[d.haplotype == h]
    pref = sub.prefix.value_counts()
    span = (sub.acc_num.max() - sub.acc_num.min()) if sub.acc_num.notna().any() else np.nan
    top_frac = pref.iloc[0] / len(sub)
    verdict = ("LIKELY SINGLE BATCH" if top_frac > 0.8 and span < 5000
               else "spread across deposits" if pref.size >= 3 else "check manually")
    batch.append({"haplotype": h, "n": len(sub), "change": haplo.set_index('haplotype')
                  .differences_vs_dominant[h], "n_prefixes": pref.size,
                  "top_prefix": pref.index[0], "top_prefix_frac": round(top_frac, 3),
                  "acc_span": span, "years": sorted(sub.year.dropna().unique().astype(int)),
                  "verdict": verdict})
batch = pd.DataFrame(batch)
batch.to_csv(OUT / "batch_diagnostic.csv", index=False)
log(batch.to_string(index=False))

# temporal denominators
tmp = []
for h in batch.haplotype:
    sub = d[d.haplotype == h]
    grp = sub.clade_imputed.mode()[0]
    for y in sorted(sub.year.dropna().unique()):
        denom = ((d.clade_imputed == grp) & (d.year == y)).sum()
        tmp.append({"haplotype": h, "year": int(y), "n": int((sub.year == y).sum()),
                    "denominator": int(denom),
                    "percent_of_clade_year": round(100 * (sub.year == y).sum() / denom, 2) if denom else np.nan})
temporal = pd.DataFrame(tmp)
temporal.to_csv(OUT / "temporal_frequencies.csv", index=False)
log("\n[13] temporal frequency (denominator = same imputed clade group, same year)")
log(temporal.to_string(index=False))

d.drop(columns=["seq", "bad_aa"]).to_csv(OUT / "per_sequence.csv", index=False)

# ============================== 7 · FIGURES ================================ #
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "axes.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})

fig = plt.figure(figsize=(7.4, 5.8))
gs  = fig.add_gridspec(2, 2, hspace=0.80, wspace=0.42)
axA = fig.add_subplot(gs[0, 0])
gsB = gs[0, 1].subgridspec(2, 1, height_ratios=[6, 1], hspace=0.05)
axB, axBd = fig.add_subplot(gsB[0]), fig.add_subplot(gsB[1])
axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])

# --- (a) clade x residue pair -------------------------------------------- #
mat = ct_pair.reindex(index=["I", "II"]).fillna(0)
axA.imshow(np.log10(mat.values + 1), cmap="Blues", aspect="auto", vmin=0)
axA.set_xticks(range(mat.shape[1])); axA.set_xticklabels(mat.columns, fontsize=8)
axA.set_yticks(range(mat.shape[0])); axA.set_yticklabels(["Clade I", "Clade II"], fontsize=8)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = int(mat.values[i, j])
        axA.text(j, i, f"{v:,}", ha="center", va="center", fontsize=11,
                 color="white" if np.log10(v + 1) > 2 else ("#B00020" if v == 0 else "0.15"),
                 fontweight="bold")
axA.text(-0.17, 1.16, "a", transform=axA.transAxes,
         fontsize=11, fontweight="bold", va="top", ha="left")
axA.text(0, 1.12, "Residues 74/107 discriminate clade", transform=axA.transAxes,
         fontsize=8.5, fontweight="bold", va="bottom")
axA.set_xlabel(f"{disc} discordant of {len(lab):,} clade-assigned sequences",
               fontsize=7, color="#B00020", labelpad=6)
for sp in axA.spines.values(): sp.set_visible(True)

# --- (b) variability along the protein, with a domain strip below --------- #
hgt = np.zeros(len(REF))
for _, r in var.iterrows(): hgt[r.position - 1] = r.minor_n
axB.bar(range(1, len(REF) + 1), hgt, width=1.0, color="#0072B2", zorder=3)
axB.set_yscale("log"); axB.set_ylim(0.7, hgt.max() * 22)
axB.set_xlim(0.5, len(REF) + .5); axB.set_xticklabels([])
axB.set_ylabel("Minor-allele count")
axB.text(-0.17, 1.16, "b", transform=axB.transAxes,
         fontsize=11, fontweight="bold", va="top", ha="left")
axB.text(0, 1.12, "Variation across 110 residues", transform=axB.transAxes,
         fontsize=8.5, fontweight="bold", va="bottom")
for _, r in var.iterrows():
    if r.minor_n >= 15:
        axB.annotate(f"{r.reference_aa}{r.position}{r.minor.split('=')[0]}",
                     xy=(r.position, r.minor_n * 1.5), ha="center", va="bottom",
                     fontsize=6.5, rotation=90)
for s0, e0, lab_, col in DOMAINS:
    axBd.add_patch(Rectangle((s0 - .5, 0), e0 - s0 + 1, 1, facecolor=col, edgecolor="white", lw=.6))
    if e0 - s0 >= 7:
        axBd.text((s0 + e0) / 2, .5, lab_, ha="center", va="center", fontsize=6, color="0.25")
axBd.set_xlim(0.5, len(REF) + .5); axBd.set_ylim(0, 1); axBd.set_yticks([])
axBd.set_xlabel("A29 residue")
for sp in axBd.spines.values(): sp.set_visible(False)

# --- (c) temporal frequency ----------------------------------------------- #
for k, hpl in enumerate(temporal.haplotype.unique()):
    t = temporal[temporal.haplotype == hpl].sort_values("year")
    chg = batch.set_index("haplotype").change[hpl]
    axC.plot(t.year, t.percent_of_clade_year, "o-", ms=4.5, lw=1.5,
             color=OKABE[k % len(OKABE)], label=f"{hpl}  {chg}")
    for _, r in t.iterrows():
        axC.annotate(f"{r.n}/{r.denominator}", xy=(r.year, r.percent_of_clade_year),
                     xytext=(0, 6), textcoords="offset points", ha="center",
                     fontsize=5.5, color=OKABE[k % len(OKABE)])
axC.set_xlabel("Year"); axC.set_ylabel("% of same-clade sequences")
axC.set_ylim(bottom=-0.8)
axC.legend(fontsize=6, frameon=False, loc="upper left", handlelength=1.2)
axC.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
axC.text(-0.17, 1.16, "c", transform=axC.transAxes,
         fontsize=11, fontweight="bold", va="top", ha="left")
axC.text(0, 1.12, "Low-frequency variants over time", transform=axC.transAxes,
         fontsize=8.5, fontweight="bold", va="bottom")

# --- (d) batch diagnostic: deposit diversity ------------------------------ #
b = batch.iloc[::-1]
y = np.arange(len(b))
axD.barh(y, b.top_prefix_frac * 100, color=[OKABE[(len(b)-1-i) % len(OKABE)] for i in range(len(b))],
         height=.55, zorder=3)
axD.axvline(80, color="#B00020", ls="--", lw=1)
axD.text(80, len(b) - .35, " single-deposit\n threshold", fontsize=6, color="#B00020", va="top")
axD.set_yticks(y); axD.set_yticklabels([f"{r.haplotype}  n={r.n}" for r in b.itertuples()], fontsize=7)
axD.set_xlim(0, 108); axD.set_xlabel("% of sequences in the largest single deposit")
for i, r in enumerate(b.itertuples()):
    axD.text(r.top_prefix_frac * 100 + 2.5, i, f"{r.n_prefixes} deposit" + ("s" if r.n_prefixes != 1 else ""), va="center", fontsize=6.5)
axD.text(-0.17, 1.16, "d", transform=axD.transAxes,
         fontsize=11, fontweight="bold", va="top", ha="left")
axD.text(0, 1.12, "Submission-batch diagnostic", transform=axD.transAxes,
         fontsize=8.5, fontweight="bold", va="bottom")

for f in ("pdf", "png", "svg"):
    fig.savefig(OUT / f"a29l_haplotype_analysis.{f}",
                dpi=600 if f == "png" else None, bbox_inches="tight")
plt.close(fig)

(OUT / "qc_report.txt").write_text("\n".join(log_lines), encoding="utf-8")
print(f"\nWrote to {OUT}/")
for f in sorted(OUT.glob("*")): print("   ", f.name)
print("\nREMINDER: haplotype frequencies measure sequencing effort, not prevalence.")
