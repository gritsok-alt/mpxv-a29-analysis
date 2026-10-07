# Analysis of the MPXV A29 predictions

Everything downstream of the predictions in `7_structure_prediction/`: topology
classification, comparison with the vaccinia A27 crystal structure, ensemble
variability, coiled-coil packing, solvent accessibility, stereochemical
validation and energy minimisation.

Three subfolders document themselves and are only summarised here; follow the
links for the detail.

**The reported model is `A29_s2m0`** — seed 2, model 0 of the reported ensemble
(`../4_reported_model/results_MPXV_A29_S1-f_seed2/`). It is the structure used by
folders 5 to 8.

**Chain convention.** Chain identifiers assigned by a predictor are arbitrary,
because the three chains are identical in sequence. Every analysis here first
relabels them so that **C is the antiparallel chain, A presents the heptad *g*
face to C, and B presents the *e* face**, matching the crystal structure. The
relabelling is determined from helix-axis dot products and Cβ–Cβ contacts at *e*
and *g* positions, never by sequence alignment. Comparisons between structures
are meaningless without it.

**Residue numbering** is native A27/A29 throughout; the 3VOP crystal is
renumbered by +20 before any comparison.

---

## 1. `1_comparison_a29_vs_a27_by_model/`

Whether the two orthologues are predicted to adopt the same structure by the
same method. `comparison_A27_A29_best_model_by_tool.pse` holds the best A29 and
A27 model from each of the three predictors; `output_RMSD.txt` records the Cα
RMSD between the two over residues 44–100, all three chains, computed in PyMOL
3.1.8 with `align ... cycles=0`.

`cycles=0` matters: PyMOL's default of five refinement cycles discards the
worst-fitting atom pairs and reports an RMSD over a subset, which is always
lower.

---

## 2. `2_chai1_coordinates/`

Supports the exclusion of Chai-1 from the benchmark.
`compare_chai1_coordinates.py` compares structures pairwise by Cα coordinates,
reporting both the maximum deviation as written and the RMSD after
superposition, which distinguishes byte-identical files from the same structure
in a different frame and from genuinely different models.
`chai1_comparison_a27.txt` and `chai1_comparison_a29.txt` are its output.

Reasoning and conclusion: [`../1_unrestrained/excluded_chai1/README.md`](../1_unrestrained/excluded_chai1/README.md).

---

## 3. `3_restraints_comparison_vs_3vop/`

How close each restraint set brings the model to the crystal arrangement.

`compare_sets_vs_3vop.py` prepares the crystal (non-polymer and zero-occupancy
atoms removed, numbering shifted +20), classifies each model's topology from its
helix axes, skips any model that is not antiparallel, relabels the rest to the
crystal convention, and reports three measurements over residues 45–84 and
47–65: assembly RMSD with all three chains fitted together, per-chain RMSD with
each chain fitted alone, and anchored deviations with one chain fitted and the
others measured in place. The three answer different questions — fold accuracy
is per-chain, packing accuracy is anchored.

```bash
python compare_sets_vs_3vop.py ../2_restraints_screening \
    --xtal ../2_restraints_screening/3vop-assembly1.cif \
    --out comparison_vs_3vop
```

It writes `comparison_vs_3vop/per_model.tsv` (one row per model per residue
window) and `comparison_vs_3vop/per_set.tsv` (the range across the models of each
set). Requires numpy and gemmi.

`pymol_session_log_s2m0_s10m0_vs_3vop.txt` is a separate, interactive comparison
of two reported-ensemble models against the crystal, run in PyMOL 3.1.8 on
24 September 2026. The script covers seed 1 of each restraint set; the log covers
these two models, which the script does not touch. The commands in the log are
reproducible by pasting them into PyMOL. `evaluations_vs_3vop-s2m0-s10m0.pse` is
the corresponding session.

---

## 4. `4_ensemble_a29_analysis/`

The step the rest of the analysis depends on: it defines the accepted ensemble
and selects the reported model.

`analyse_ensemble.py` classifies chain orientation, relabels chains to the common
convention, checks each imposed restraint against its 8 Å bound, determines the
heptad face each parallel chain presents, computes the per-residue spread across
accepted models, and selects the reported model. Acceptance requires the
two-parallel/one-antiparallel arrangement and ipTM ≥ 0.60; **97 of the 100 models
were accepted**. Residue numbering, restraint set, distance bound and helix-fit
range are command-line options, with examples in the docstring. Requires numpy
and gemmi.

Outputs in `a29_reported/`:

| File | Contents |
| --- | --- |
| `per_model.tsv` | one row per model: topology class, chain map, restraint distances, confidence scores |
| `per_seed.tsv` | one row per seed |
| `per_residue.tsv` | per-residue RMSF and mean pLDDT across the accepted models |
| `analysis.xlsx` | the same three tables as worksheets |
| `glossary.txt` | the column definitions |
| `relabelled/` | the accepted models after chain relabelling |

`RMSF_plot.xlsx` is the RMSF figure, drawn manually from `per_residue.tsv`. The
spreadsheet is a convenience copy; the TSV is the source.

---

## 5. `5_socket2/`

Independent validation of the heptad register: Socket v3.02 assigns the register
from coordinates by locating knobs-into-holes packing, and the result is compared
with the register derived from sequence before any model existed. Both structures
and both packing cutoffs give an uninterrupted register in agreement.

[`5_socket2/README.md`](5_socket2/README.md).

---

## 6. `6_solvent_accessible_surface_area/`

Accessibility of functionally annotated residues across the accepted ensemble,
testing the model against published mutagenesis that was never supplied to the
prediction. Reported as distributions over 97 models rather than from a single
structure, with an occlusion control in which residues 1–43 are removed.

[`6_solvent_accessible_surface_area/README.md`](6_solvent_accessible_surface_area/README.md).

---

## 7. `7_molprobity/`

Stereochemical validation on the MolProbity web server, one folder per input:
`3vop/` (the crystal structure), `A29_s2m0/` (the reported model as predicted),
and `A29_s2m0_minimized_openmm/` and `A29_s2m0_minimized_yasara/` (after the two
minimisations in folder 8).

Every input was treated identically: existing hydrogens removed, then added by
MolProbity with Asn/Gln/His flips at electron-cloud rather than nuclear
positions — the `FH` in each output filename. Uniform treatment is what makes the
clashscores comparable across the four. `sequence.docx` in each folder records
the steps taken.

**The validation measurements reported in the article — burial, register, RMSD to
the crystal — were made on the unminimised coordinates.** Minimisation was
applied only so that a relaxed structure could be deposited.

---

## 8. `8_energy_minimization/`

Two independent minimisations of the reported model, for comparison.

`openmm/minimise.py` runs OpenMM 8.6.1 with the Amber14 force field in vacuo,
after PDBFixer preparation (hydrogens rebuilt at pH 7.4, unresolved loops not
rebuilt, heterogens and water removed). Backbone N, CA, C and O are restrained at
10 kcal/(mol·Å²) so that side chains relax without the fold moving, and
minimisation runs to a tolerance of 10 kJ/(mol·nm). The script reports the RMSD
from the input coordinates, which is the evidence that the fold did not shift.
Needs its own environment (OpenMM 8.6.1, PDBFixer 1.12).

`yasara/` holds the result of the YASARA energy-minimisation web server, as the
scene file and as the PDB exported from it with YASARA View 26.6.12.W.64.

---

## Environments

| Script | Needs |
| --- | --- |
| `compare_sets_vs_3vop.py`, `analyse_ensemble.py` | numpy, gemmi |
| `ensemble_sasa.py` | numpy, gemmi, Biopython |
| `compare_chai1_coordinates.py` | numpy, Biopython |
| `minimise.py` | OpenMM, PDBFixer |

The first three were run in the `r4s` conda environment (numpy 2.5.3,
gemmi 0.7.5, Biopython 1.88); `minimise.py` in a separate `openmm` environment
(OpenMM 8.6.1, PDBFixer 1.12).

---

## Helix-axis fitting ranges

Three appear in this project and are not interchangeable. `analyse_ensemble.py`
and `ensemble_sasa.py` fit residues 44–100; `compare_sets_vs_3vop.py` and the
PyMOL session log fit 47–65; the prediction notebooks' in-session topology check
uses 44–84. All three classify orientation, and a borderline model could in
principle classify differently under each. **The 44–100 range is the one the
reported acceptance set and RMSF rest on.**


