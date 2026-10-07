# Step 6 — Secondary architecture of MPXV A29

This folder contains the sequence-based predictions of the secondary
architecture of MPXV A29: signal peptide and transmembrane topology, secondary
structure and intrinsic disorder, coiled-coil propensity and heptad register,
and the physicochemical properties of the long helix.

All analyses were run on public **web services**, so there is no code in this
folder; the downloaded result files are the permanent record. Every prediction
used the MPXV A29L reference sequence **YFT67316.1** (110 aa), deposited here as
`Q77HM6.fasta`. Some result files label the query with the UniProt accession
**Q77HM6**, whose sequence is identical.

Residue numbering throughout follows YFT67316.1. Domain boundaries used for
interpretation, matching Tables S4 and S5: NTR 1–20, GAG-binding domain 21–32,
spacer 33–42, CCD 43–84, LZD 85–110.

---

## Folder layout

```
6_secondary_architecture/
├── README.md
├── Q77HM6.fasta            the query sequence (YFT67316.1 = Q77HM6)
├── A29L_msa_linsi.fasta    60-sequence alignment, for the PCOILS MSA runs
├── topology/
│   ├── signalp/            SignalP 6.0
│   ├── phobius/            Phobius
│   ├── deeptmhmm/          DeepTMHMM
│   └── memsat_svm/         MEMSAT-SVM (PSIPRED Workbench)
├── structure_disorder/
│   ├── psipred/            PSIPRED 4.0 (PSIPRED Workbench)
│   ├── s4pred/             S4PRED (PSIPRED Workbench)
│   ├── disopred/           DISOPRED3 (PSIPRED Workbench)
│   └── iupred3/            IUPred3 and ANCHOR2
├── coiled_coil/
│   ├── marcoil/            MARCOIL (MPI Bioinformatics Toolkit)
│   ├── deepcoil2/          DeepCoil2 (MPI Bioinformatics Toolkit)
│   ├── coconat/            CoCoNat
│   └── pcoils/
│       ├── single_sequence/   PCOILS on the reference sequence
│       └── msa/               PCOILS on the 60-sequence alignment
└── helical_properties/
    ├── heliquest/          HeliQuest
    └── waggawagga/         WaggaWagga
```

`A29L_msa_linsi.fasta` is a copy of
`../1_orthologs_msa/output_phase7/A29L_msa_linsi.fasta`, kept here because the
PCOILS MSA runs cannot be reproduced without it.

---

## 1. Signal peptide and transmembrane topology (`topology/`)

| Tool | Service | Settings | Files dated |
| --- | --- | --- | --- |
| SignalP 6.0 | DTU web server | Eukarya; long output; slow model | 2 Sep 2026 |
| Phobius | EMBL-EBI Job Dispatcher | long-run mode, graphics output | 22 Sep 2026 |
| DeepTMHMM | BioLib web server | default | 22 Sep 2026 |
| MEMSAT-SVM | PSIPRED Workbench | default | 31 Jul 2026 (job `fe6a07da…`) |

The predictions were first made on 31 July 2026. SignalP, Phobius and DeepTMHMM
were later rerun, giving identical predictions; the saved files are from those
reruns, and the 31 July outputs were **not** retained, so the identity of the two
sets cannot be verified from this deposit.

**Results.** No signal peptide was detected by any method: SignalP classified the
sequence as OTHER (probability 1.00, 0 for every signal-peptide class),
DeepTMHMM assigned it to the globular class (GLOB), and Phobius and MEMSAT-SVM
likewise detected none. The empty `signalp_processed_entries.fasta` and the
header-only GFF3 files are the direct expression of this: SignalP found no
cleavage site to report. These predictions do not support the historical
annotation of residues 1–20 as a signal peptide, inherited from the vaccinia A27
literature.

On membrane insertion, DeepTMHMM and Phobius predicted no transmembrane segment,
whereas MEMSAT-SVM predicted a single helix at residues 91–106, classified as
pore-lining (posterior probability 0.986) with a predicted pore stoichiometry of
four. This assignment is weakly supported: every multi-helix topology (two to
five helices, both orientations) was rejected during the topology search, the two
orientations of the single helix scored almost identically (−0.7381 and −0.7361),
and the segment contains charged residues incompatible with a membrane-spanning
helix. Residues 91–106 coincide with the hydrophobic core of the coiled coil in
the LZD (see section 3), which accounts for its periodic hydrophobicity.

**Files.** `signalp/`: `signalp_prediction_results.txt` (summary),
`signalp_output.json`, `signalp_output_sequence_plot.*` (per-residue
probabilities as text, PNG and EPS) and GFF3 region files;
`signalp_processed_entries.fasta` is empty because no signal peptide was cleaved.
`phobius/`: `phobius.out`, `phobius_plot.png`. `deeptmhmm/`:
`deeptmhmm_predicted_topologies.3line`, `TMRs.gff3`, `plot_deeptmhmm.png`.
`memsat_svm/`: see Figure S9 below.

**Supplementary Figure S9** is assembled from the MEMSAT-SVM job:
panel A = `memsat_schematic.png`, B = `cartoon_memsat_svm.png`,
C = `annotation_grid_2.png` (MEMSAT-SVM prediction mapped on the sequence),
D = `annotation_grid_1.png` (residue classes only; not a prediction).
The complete MEMSAT-SVM output, including all rejected topologies, is
`memsat_results.memsat_svm`.

---

## 2. Secondary structure and intrinsic disorder (`structure_disorder/`)

| Tool | Service | Settings | Files dated |
| --- | --- | --- | --- |
| PSIPRED 4.0 | PSIPRED Workbench | default | 20–21 Aug 2026 |
| S4PRED | PSIPRED Workbench | default | 21 Aug 2026 |
| DISOPRED3 | PSIPRED Workbench | default | 20 Aug 2026 |
| IUPred3, ANCHOR2 | IUPred3 web server | short and long modes | 31 Jul 2026 |

**Files.** `psipred/`: `psipred.ss2`, `psipred.horiz`, `psipred_chart.svg`,
`psipred_annotation_grid.svg`. `s4pred/`: `s4pred.ss2`, `s4pred.horiz`.
`disopred/`: `disopred.comb` (per-residue disorder), `disopred.pbdat`
(protein-binding disorder), `disopred_chart.png`. `iupred3/`:
`iupred3_short_disorder.txt`, `iupred3_long_disorder.txt`, `anchor2.txt`.

**Secondary structure.** PSIPRED and S4PRED gave residue-for-residue identical
predictions: coil 1–28, a short helix 29–38, coil 39–43, a long α-helix 44–100,
and coil 101–110 (`.horiz` files give the confidence per residue).

**Disorder.** Disorder is predicted sparingly and largely at the termini.
DISOPRED3 assigns disorder to residues 1–4 only; its mean score falls from 0.424
across the NTR to 0.041 across the CCD. IUPred3 in short mode exceeds 0.5 only at
residues 1–9 and 104–110. In long mode, which is parameterised for disordered
regions of at least 30 residues, no run above 0.5 exceeds nine residues (the
longest being 32–40, spanning the end of the GAG-binding domain and the spacer).
ANCHOR2 stays below 0.5 throughout (maximum 0.486), so no disordered binding
region is predicted. **The article reports IUPred3 in short mode**;
`iupred3_long_disorder.txt` supports the statement on long mode, and
`anchor2.txt` contains both the IUPred3 long-mode and ANCHOR2 scores.

---

## 3. Coiled coil and heptad register (`coiled_coil/`)

| Tool | Service | Settings | Files dated |
| --- | --- | --- | --- |
| MARCOIL | MPI Bioinformatics Toolkit | default | 13 Aug 2026 |
| DeepCoil2 | MPI Bioinformatics Toolkit | default | 13 Aug 2026 |
| CoCoNat | CoCoNat web server | default | 14 Aug 2026 |
| PCOILS | MPI Bioinformatics Toolkit | MTIDK matrix, sequence weighting, windows 14/21/28 | 14 Aug 2026 |

PCOILS was run on the single reference sequence (`pcoils/single_sequence/`) and
on the 60-sequence orthologue alignment (`A29L_msa_linsi.fasta`; `pcoils/msa/`),
both with sequence weighting. Each run reports three window sizes, giving six
PCOILS configurations.

**Files.** `marcoil/marcoil.out` (per-residue probability and heptad phase);
`deepcoil2/deepcoil2.out` and `img_deepcoil2.png`;
`coconat/results_coconat.txt` (segments, register and oligomeric state),
`coconat.tsv` (per-residue probabilities) and `coconat.png`;
`pcoils/single_sequence/ss-pcoils.out` and `pcoils/msa/msa-pcoils.out`
(per-residue probability and register for all three windows), each with its
coiled-coil (`*img_ncoils.png`) and secondary-structure (`*img_psipred.png`)
plots.

**Extent of the coiled coil.** All methods place a single coiled coil within the
long helix, and none in the N-terminal region, but the boundaries differ:

| Method | Coiled-coil region | Notes |
| --- | --- | --- |
| MARCOIL | 40–74 (probability > 90%) | drops to 16% at residue 75, below 5% through the LZD |
| DeepCoil2 | 46–97 (probability > 0.5; plateau 0.78) | continuous |
| CoCoNat | 45–101 | continuous; predicted oligomeric state tetramer (probability 0.65) |
| PCOILS, single sequence, windows 14 / 21 / 28 | 42–100 / 41–102 / 41–108 | each window interrupted, at residues 64 / 70 / 74 |
| PCOILS, alignment, windows 14 / 21 / 28 | 48–101 / 41–103 / 41–109 | only window 14 interrupted (residue 64); windows 21 and 28 continuous |

Including homologues removes the interruption at the longer windows: the
single-sequence runs break at residues 70 and 74 for windows 21 and 28, whereas
the corresponding alignment runs are continuous. The window-14 break at residue
64 persists in both. The interruptions fall across the polar core pair at the
72–78 heptad (Cys72 at *a*, Asn75 at *d*), the same region where MARCOIL
terminates — so the disagreement between methods is concentrated at one
structurally awkward position rather than distributed along the helix.

Cys72 is also one half of the disulfide-linked motif at residues 71–72, which
Tables S4 and S5 report as invariant across both the 60-orthologue and the
14-orthologue sets. A position that the coiled-coil predictors find difficult is
therefore the most strongly conserved in the alignment.

**Heptad register.** The register is identical across MARCOIL, CoCoNat, DeepCoil2
and all six PCOILS configurations: *a* = residues 51, 58, 65, 72, 79, 86 and 93;
*d* = residues 47, 54, 61, 68, 75, 82, 89 and 96.

The article's coiled-coil confidence figure overlays the per-residue profiles
from MARCOIL, CoCoNat, DeepCoil2 and PCOILS window 21 (single sequence and
alignment), all taken from the files listed above.

---

## 4. Helical properties (`helical_properties/`)

**HeliQuest** (14 Aug 2026; α-helix, 18-residue sliding window, 93 windows).
`heliquest.txt` lists, for every window, the mean hydrophobicity ⟨H⟩, the
hydrophobic moment ⟨µH⟩, the net charge and the polar/non-polar fractions;
`plot_mom_hyd.jpeg` plots ⟨H⟩ and ⟨µH⟩ along the sequence. The hydrophobic
moment peaks within the long helix, at 0.567 for window 51–68 and 0.553 for
window 65–82, where HeliQuest identifies a continuous hydrophobic face
(e.g. I-I-L-F-I-L). Windows beginning before residue 39 never exceed 0.32.

`heliquest.pdf` (16 MB) is a print of the HeliQuest result page, showing every
window with its helical wheel. The service offers no other export of the wheels,
so the page print is the only available record; its size is due to the embedded
per-window images.

**WaggaWagga** was used as a visual guide for the heptad organisation; the heptad
figure in the article was redrawn manually in Inkscape. The saved views
(`helical_net_view.png`, `helical_wheel_view.png`, `SAH.svg`) come from
WaggaWagga's internal MARCOIL prediction, which spans residues 51–99 and shows
the same register as above. Note that this segment differs from the standalone
MARCOIL run in `coiled_coil/marcoil/` (40–74), since WaggaWagga runs its own
implementation and settings. WaggaWagga's single α-helix (SAH) score for this
segment is 0.049, i.e. the helix is not predicted to be a stable single α-helix,
consistent with a coiled coil.
