# Step 2 — Maximum-likelihood phylogeny of A29L

This folder contains the phylogenetic analysis described in section 2.1.3 of the
article: a maximum-likelihood gene tree of the A29L family, branch support,
per-genus monophyly tests, and a substitution-saturation diagnostic.

**This is a gene tree, not a species tree.** It describes the evolutionary history
of the A29L protein, which need not match the history of the viruses carrying it.
A29L is only ~110 aa, so deep nodes are expected to be weakly supported — and the
saturation test below measures why.

---

## Folder layout

```
2_phylogeny/
├── README.md          this file
├── IQTREE.ipynb       tree search, branch support, monophyly tests
├── fig_S7.png         iTOL tree after improvement in InkScape (Fig. S7) 
├── output_iqtree/     all IQ-TREE results, tables, logs and the tree figure
└── saturation/        substitution-saturation diagnostic (Fig. S8)
```

## Input

The notebook reads the final alignment produced in Step 1 directly:

```
../1_orthologs_msa/output_phase7/A29L_msa_linsi.fasta   (60 sequences × 300 columns)
```

It also imports the shared FASTA parser `seqio_utils.py` from
`../1_orthologs_msa/`. No copy of either file is kept beside the notebook, so
Step 1 must be present beside this folder for it to run. (`saturation/` keeps its
own copies of both inputs; see that folder's README.)

Sequence headers follow the `Genus|Accession|Species` convention from Step 1, so
genus labels for colouring and for the monophyly tests are read directly from the
alignment.

---

## Software

| Software | Version | Notes |
|---|---|---|
| IQ-TREE | 3.0.1 (Linux x86 64-bit) | called as `iqtree3` through WSL |
| Python | 3.12 | `matplotlib`; `numpy` for the saturation test |
| iTOL | 7.6 | final rooting and visualisation (web service) |
| Inkscape | 1.3.2 | final figure refinement |
| Microsoft Excel | — | Fig. S8 drawn manually; see `saturation/` |

On a native Linux system, set `wsl_prefix` to `[]` in the notebook's `Config`.

---

## Analyses

Run the notebook from top to bottom with this folder as the working directory.

**1. Model selection, tree search and branch support** in a single IQ-TREE run:

```
iqtree3 -s A29L_msa_linsi.fasta -m MFP -B 1000 --alrt 1000 --abayes -T AUTO --seed 12345 --prefix output_iqtree/A29L
```

ModelFinder (`-m MFP`) selects the best-fit model by BIC. Branch support is
estimated by 1,000 ultrafast bootstrap replicates (UFBoot), 1,000 SH-aLRT
replicates and the approximate Bayes test (aBayes). Support values appear on each
internal node in the order **SH-aLRT / aBayes / UFBoot**. A branch is considered
well supported when **UFBoot ≥ 95 and SH-aLRT ≥ 80**.

`-T AUTO` selected **4 threads** on the machine used. The seed reproduces this run
exactly only at the same thread count: IQ-TREE consumes random numbers in a
different order under a different degree of parallelism, so a reader on other
hardware may obtain a slightly different tree despite `--seed 12345`.

**2. Per-genus monophyly tests.** For each genus with at least two sequences, a
constraint tree forcing only that genus to be monophyletic is searched under the
best-fit model (`-g`), and the constrained tree is compared with the unconstrained
ML tree by the approximately unbiased (AU) test (`-z … -n 0 -zb 10000 -au`).
Sequences labelled `unclassified` are excluded, as they share a missing label
rather than a genus; genera with a single sequence are trivially monophyletic and
are not tested.

**3. Tree figure.** A genus-coloured tree, midpoint-rooted for display only, with
well-supported branches marked (`A29L_ml_tree.png`, 1000 dpi).

The whole run is checkpointed: if `output_iqtree/A29L.treefile` exists, the tree
search is skipped, and each AU test is skipped if its report exists. To reproduce
the analysis from scratch, delete `output_iqtree/`.

---

## Substitution saturation

`saturation/` holds a diagnostic applied to the tree and alignment above: pairwise
observed p-distance plotted against pairwise patristic distance, against the
p-distance expected between sequences of the same amino-acid composition (0.93).
Beyond about 2 substitutions per site the relationship decelerates sharply — the
binned means rise only from 0.58 to 0.77 across the remaining five-fold increase
in patristic distance — and 76% of all taxon pairs lie in that region.

This is the mechanism behind the weak backbone support reported below and behind
the AU tests' inability to distinguish constrained from unconstrained topologies:
most pairs in the dataset carry little residual information about deep branching
order. See `saturation/README.md`. The figure (Fig. S8) was drawn in Excel from
the script's CSV output; the script itself produces statistics, not a figure.

---

## Results

| Item | Value |
|---|---|
| Alignment | 60 sequences × 300 columns, 14 genera + 4 unclassified sequences |
| Best-fit model (BIC) | **Q.INSECT+F+I+G4** |
| BIC score | 14808.5149 |
| Internal branches | 57 |
| Branches with UFBoot ≥ 95 | 18 |
| Branches with SH-aLRT ≥ 80 | 26 |
| Well-supported branches (both criteria) | 16 (28.1%) |

Per-genus monophyly (AU test against the unconstrained ML tree):

| Genus | Sequences | ΔlogL | p-AU | Consistent with monophyly |
| --- | --- | --- | --- | --- |
| Centapoxvirus | 3 | 3.21 | 0.382 | yes |
| Leporipoxvirus | 4 | 2.70 | 0.383 | yes |
| Yatapoxvirus | 2 | 2.70 | 0.390 | yes |
| Oryzopoxvirus | 2 | 1.25 | 0.402 | yes |
| Orthopoxvirus | 14 | 1.24 | 0.405 | yes |
| Parapoxvirus | 16 | 0.92 | 0.427 | yes |

All six tested genera are consistent with monophyly (p-AU > 0.05). The p-AU values
fall in a narrow band (0.38–0.43) and none approaches rejection, which on a ~110 aa
alignment reflects limited power as much as positive agreement: non-rejection here
means the data cannot distinguish the constrained from the unconstrained topology,
not that monophyly is well supported. The saturation test above quantifies that
limitation directly.

Cervidpoxvirus (4 sequences) and Capripoxvirus (5) are monophyletic in the
unconstrained ML tree, so the constrained search recovered the ML topology itself
(ΔlogL = 0.00) and no comparison was possible. They are reported as monophyletic
without an AU test; the p-AU values their runs returned compare a tree with
itself and are not interpretable.

Not tested: `unclassified` (4 sequences), and the single-sequence genera
Avipoxvirus, Mustelpoxvirus, Pteropopoxvirus, Sciuripoxvirus, Suipoxvirus and
Vespertilionpoxvirus.

This tests whether the history of the A29L protein is consistent with the genus
assignments; it is not a test of poxvirus taxonomy.

---

## Rooting and the article figure

IQ-TREE produces unrooted trees, and the notebook's own figure is midpoint-rooted
for display only. The tree shown in the article was produced from
`output_iqtree/A29L.treefile` in **iTOL v7.6**, rooted on the single Avipoxvirus
sequence, and refined in **Inkscape v1.3.2**.

---

## Output files

| File | Content |
| --- | --- |
| `A29L.treefile` | ML tree with support values (Newick); basis of the article figure and of Step 3 |
| `A29L.contree` | UFBoot consensus tree |
| `A29L.iqtree` | Full report: model ranking, parameters, tree statistics |
| `A29L.log` | IQ-TREE log of the main run |
| `A29L.model.gz`, `A29L.mldist`, `A29L.bionj`, `A29L.splits.nex` | Model-selection results, ML distances, starting tree, split supports |
| `branch_support.csv` | Support values for every internal branch |
| `genus_monophyly_tests.csv` | AU-test results per genus |
| `constraint_<genus>.nwk` | Constraint used for each genus |
| `constrained_<genus>.*` | Best tree under each constraint |
| `candidates_<genus>.nwk` | ML and constrained trees compared in each AU test |
| `AU_<genus>.*` | AU-test reports |
| `A29L_ml_tree.png` | Genus-coloured tree figure (midpoint-rooted) |
| `iqtree_analysis.log` | Log of the notebook itself |

Files named `iqtree_run.log` and `autest_<genus>.log` are the notebook's captures
of IQ-TREE's screen output and duplicate the corresponding IQ-TREE logs. IQ-TREE
checkpoint files (`*.ckp.gz`) are not included.

---

## Note on labels

Tip labels are taken from the Step 1 alignment headers. The label
`Orthopoxvirus_ectomelia` is a misspelling of *ectromelia* (ectromelia virus)
inherited from the Step 1 curation. It does not affect any result, and was
corrected by hand in iTOL for the article figure — so the published figure and the
deposited tree file differ in this one label.
