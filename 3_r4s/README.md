# Site-specific evolutionary rates (Rate4Site)

Per-residue evolutionary rates for MPXV A29 (OPG154), referenced to YFT67316.1
(110 aa), inferred at two taxonomic levels:

| level | alignment | sequences | fitted α |
| --- | --- | --- | --- |
| Chordopoxvirinae | `60sequences_subfamily.fasta` | 60 | 1.35 |
| Orthopoxvirus | `14sequences_genus.fasta` | 14 | 0.834 |

Two further runs at the Chordopoxvirinae level, under alternative substitution
matrices, support the robustness check reported in the manuscript.

**Date run:** 1 October 2026 (JTT runs); 6 October 2026 (WAG, LG)

---

## Method

Rate4Site infers the relative evolutionary rate at each alignment column by
empirical Bayes, conditioned on a fixed tree topology, under a gamma
distribution of rates across sites. The output is normalised to mean 0 and
standard deviation 1 across the reference positions, so that **lower scores
indicate stronger conservation**.

| | |
| --- | --- |
| Program | Rate4Site v3.0.0 |
| Substitution matrix | JTT (`-Mj`); WAG (`-Mw`) and LG (`-Ml`) for the robustness check |
| Rate inference | empirical Bayes (`-ib`) |
| Gamma categories | 16 (`-k 16`) |
| Reference sequence | YFT67316.1 (MPXV A29, clade I) |
| Tree, Chordopoxvirinae | IQ-TREE maximum likelihood, model Q.INSECT+F+I+G4 |
| Tree, Orthopoxvirus | the same tree, pruned to the 14 *Orthopoxvirus* taxa |

Two substitution models are involved at different stages: the tree was built
under Q.INSECT+F+I+G4 as selected by IQ-TREE's ModelFinder, while Rate4Site
inferred the per-site rates under JTT, which is among the matrices it offers.
The Q-matrix series is not implemented in Rate4Site v3.0.0, so Q.INSECT rates do
not exist and are not part of any comparison here. Branch lengths were
re-optimised by Rate4Site under the gamma model; the topology was fixed to the
IQ-TREE result and never searched.

**The Orthopoxvirus topology is pruned, not re-inferred.** `prune_tree_to_subset.py`
removes the non-*Orthopoxvirus* tips from the 60-taxon maximum-likelihood tree,
adding branch lengths together when a degree-two node collapses, so distances
among the retained taxa are unchanged. Re-inferring a 14-taxon tree from a 110 aa
alignment would give a different and much less well determined topology, and the
two levels would no longer be nested. This matters: the fitted α depends on the
topology supplied, and a re-inferred tree gives a different value.

The 14-sequence alignment is a column subset of the 60-sequence alignment with
all-gap columns removed (144 columns), not a fresh MAFFT run, so column
numbering and residue correspondence are identical at both levels.

**Cite:** Pupko *et al.* (2002) *Bioinformatics* 18:S71–S77; Mayrose *et al.*
(2004) *Mol. Biol. Evol.* 21:1781–1791.

---

## Input preprocessing

Rate4Site's Newick parser fails on two features of IQ-TREE output, in both
cases silently. `run_rate4site.py` handles them:

1. **Pipe characters in taxon labels.** A label such as
   `Orthopoxvirus|YFT67316.1|Orthopoxvirus_monkeypox` is interpreted as a shell
   pipeline. All taxa are renamed to `T01`–`Tnn`; `namemap.json` restores the
   original names.
2. **Internal node labels.** IQ-TREE writes support values as
   `75.1/0.856/94` at internal nodes, which Rate4Site cannot parse. These are
   stripped. Branch lengths and topology are preserved unchanged.

The alignment and tree taxon sets are verified identical before each run
(60/60 and 14/14).

---

## Positional convention

All positions are 1–110 of the MPXV A29 reference YFT67316.1, with columns
gapped in that sequence removed (Rate4Site's
`removeUnknownPositionsAccordingToAReferenceSeq`). Both levels collapse to the
same 110 reference positions, which is what makes them directly comparable.

---

## Conservation grades and the low-confidence flag

The continuous rate is binned into nine grades following the ConSurf
convention, 9 being most conserved. Bins are **equal-interval** over the
observed score range; note that the ConSurf web server uses equal-occupancy
bins, so grades from the two sources are not directly comparable.

A position is flagged as low-confidence when its 25–75 % posterior interval
spans more than a third of the total score range — the criterion ConSurf uses
for its "insufficient data" category.

In the Chordopoxvirinae JTT run, seventeen positions were flagged: 16, 19, 21,
24, 27, 28, 30, 31, 32, 34, 35, 37, 39, 40, 91, 95 and 110. Fourteen lie between
residues 16 and 40, where alignment occupancy is sparsest; 91 and 95 are fully
occupied, so their wide intervals reflect genuine rate uncertainty rather than
missing data. These positions are plotted but should not be interpreted.

**The flagged set is approximate.** The count depends on the substitution
matrix: 17 under JTT, 18 under WAG, 22 under LG. The rate *ranking* is stable
across matrices (below), but the precision of individual estimates is not, so the
flagged positions should be treated as a region of low confidence rather than a
fixed list. The Orthopoxvirus run flags 33 of 110 positions, as expected at that
sampling depth.

---

## Reproducing

```bash
# environment
conda create -n r4s -c conda-forge biopython matplotlib pandas fonttools
conda activate r4s
# Rate4Site is installed system-wide; build from
# https://github.com/barakav/r4s_for_collab if absent
```

### Chordopoxvirinae level (60 sequences, JTT) — primary run

```bash
python run_rate4site.py \
  --msa  60sequences_subfamily.fasta \
  --tree A29L.treefile \
  --ref  YFT67316.1 \
  --out  r4s_run_subfamily_chordopoxvirus \
  --tag  r60_jtt
```

Reports `fitted gamma alpha = 1.35`.

### Orthopoxvirus level (14 sequences, JTT)

Prune the topology first, then run:

```bash
python prune_tree_to_subset.py \
  --tree A29L.treefile \
  --msa  14sequences_genus.fasta \
  --out  A29L_orthopoxvirus_pruned.treefile

python run_rate4site.py \
  --msa  14sequences_genus.fasta \
  --tree A29L_orthopoxvirus_pruned.treefile \
  --ref  YFT67316.1 \
  --out  r4s_run_orthopoxvirus \
  --tag  r14_jtt
```

Reports `fitted gamma alpha = 0.834346`.

### Substitution-matrix robustness (WAG, LG)

```bash
python run_rate4site.py \
  --msa 60sequences_subfamily.fasta --tree A29L.treefile --ref YFT67316.1 \
  --out r4s_run_subfamily_wag --tag r60_wag --matrix=-Mw

python run_rate4site.py \
  --msa 60sequences_subfamily.fasta --tree A29L.treefile --ref YFT67316.1 \
  --out r4s_run_subfamily_lg  --tag r60_lg  --matrix=-Ml
```

**Note the `--matrix=-Mw` form.** Written as `--matrix -Mw`, argparse reads the
leading hyphen as a new option and the run fails with
`expected one argument`. The equals form is required.

### Comparisons

```bash
# matrix robustness
python3 compare_rate_levels.py \
  --a r4s_run_subfamily_chordopoxvirus/r60_jtt.res --label-a JTT \
  --b r4s_run_subfamily_wag/r60_wag.res --label-b WAG

python3 compare_rate_levels.py \
  --a r4s_run_subfamily_chordopoxvirus/r60_jtt.res --label-a JTT \
  --b r4s_run_subfamily_lg/r60_lg.res --label-b LG

# cross-level
python3 compare_rate_levels.py \
  --a r4s_run_subfamily_chordopoxvirus/r60_jtt.res --label-a "60 seq (Chordopoxvirinae)" \
  --b r4s_run_orthopoxvirus/r14_jtt.res            --label-b "14 seq (Orthopoxvirus)" \
  --by-domain --csv rate_comparison.csv
```

`compare_rate_levels.py` uses only the Python standard library.

### Figure

```bash
python plot_rate4site_A29.py --res r4s_run_subfamily_chordopoxvirus/r60_jtt.res --out figure3
```

`plot_rate4site_A29.py` downloads Tinos (Apache-2.0, metrically compatible with
Times New Roman) and renames its internal records so matplotlib resolves the
family. Pass `--no-font-fetch` to use a locally installed Times instead, which is
also the offline path. The downloaded fonts are not tracked in this repository.

---

## Results of the comparisons

### Substitution matrix

Rate4Site normalises within each run, so only the ordering of positions is
comparable between runs. The ranking is essentially unaffected by the matrix:

| comparison | Spearman ρ |
| --- | --- |
| JTT vs WAG | 0.990 |
| JTT vs LG | 0.988 |

### Taxonomic level

Across all 110 positions, ρ = 0.634 (p = 1 × 10⁻¹³): the two levels agree
substantially on which residues are relatively constrained, but about 40 % of the
rank variance is not shared, so the Orthopoxvirus-level rates carry information
additional to a rescaling of the subfamily-level ones. Per domain:

| domain | residues | n | ρ | p |
| --- | --- | --- | --- | --- |
| N-terminal region | 1–20 | 20 | 0.120 | 0.614 |
| GAG-binding domain | 21–32 | 12 | 0.727 | 0.0074 |
| Spacer | 33–42 | 10 | 0.467 | 0.174 |
| Coiled-coil domain | 43–84 | 42 | 0.487 | 0.0011 |
| Leucine zipper domain | 85–110 | 26 | 0.745 | 1.3 × 10⁻⁵ |

The disulfide-linked motif (71–72) has too few positions to correlate.

**Sensitivity to the low-confidence flag.** Excluding every position flagged in
either run leaves 74 of 110 and lowers the global correlation to ρ = 0.284
(p = 0.014). This attenuation is expected rather than contradictory: wide
posterior intervals concentrate at fast-evolving positions, which are also where
the two runs most obviously agree, so removing them leaves the conserved
majority, where scores are compressed and the ordering is dominated by noise. The
agreement between levels is therefore carried substantially by the fast-evolving
positions. The unfiltered statistics are reported as primary, since filtering on a
precision criterion biases the retained set toward low variance.

A paired Wilcoxon signed-rank test on the same data is reported by the script but
is not informative here: both runs are z-scored, so a systematic shift between
them is close to excluded by construction (p = 0.48 cross-level, p = 0.79 and 0.87
for the matrix comparisons). It is retained in the output as a check, not as a
result.

---

## Files

| | |
| --- | --- |
| `60sequences_subfamily.fasta` | Chordopoxvirinae alignment (60 seq, 300 columns), from `1_orthologs_msa/output_phase7/A29L_msa_linsi.fasta` |
| `14sequences_genus.fasta` | *Orthopoxvirus* column subset (14 seq, 144 columns) |
| `A29L.treefile` | IQ-TREE maximum-likelihood tree, from `2_phylogeny/output_iqtree/` |
| `A29L_orthopoxvirus_pruned.treefile` | the same tree pruned to the 14 *Orthopoxvirus* taxa |
| `run_rate4site.py` | Preprocessing and Rate4Site invocation |
| `prune_tree_to_subset.py` | Tree pruning for the Orthopoxvirus level |
| `compare_rate_levels.py` | Spearman and Wilcoxon comparison of two `.res` files |
| `plot_rate4site_A29.py` | Figure generation |
| `r4s_run_subfamily_chordopoxvirus/` | **Primary output** — 60-sequence JTT run |
| `r4s_run_orthopoxvirus/` | 14-sequence JTT run |
| `r4s_run_subfamily_wag/`, `r4s_run_subfamily_lg/` | matrix-robustness runs |
| `rate_comparison.csv` | Per-position scores and ranks at both levels |
| `pymol_mapping_code.txt` | PyMOL commands mapping rates onto the structure |
| `figure3.png`, `figure3.pdf` | Figure 3 of the manuscript (raw) |

Each `r4s_run_*/` directory holds `<tag>.res` (rates in reference coordinates),
`<tag>.unnorm` (un-normalised rates), `<tag>.tree` (re-optimised branch lengths),
`<tag>.log`, and the sanitised inputs `aln.fasta`, `tree.nwk` and `namemap.json`.

---

## Notes

Rate4Site reports *lower* scores for more conserved positions. Figures in the
manuscript plot the score directly, so conserved positions fall below zero.

The alignment-wide mean rate is zero by construction, not an empirical quantity:
statements of the form "below the alignment mean" mean "below zero".
