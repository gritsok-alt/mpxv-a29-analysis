# Step 5 — A29 variation within circulating monkeypox virus

This folder contains the MPXV population survey described in section 2.1.6 of
the article: the variation of the A29 protein across all monkeypox virus A29
sequences available in NCBI Virus.

---

## Folder layout

```
5_mpxv_survey/
├── README.md
├── a29l_haplotypes_web.py     the complete analysis
├── sequences.fasta            NCBI Virus export (input, 9,870 sequences)
├── haplotypes17.fa            the 17 haplotype sequences (prepared in Jalview)
└── results_web/               all outputs of the script
```

## Data

`sequences.fasta` was downloaded from **NCBI Virus** on **28 July 2026**
(*Orthopoxvirus monkeypox*, A29 protein sequences of 110 aa). Headers follow
the NCBI Virus export format, including the clade and collection date:

```
>YDZ51540.1 |IMV surface fusion protein [Monkeypox virus]|IIb|2026-03-24
```

NCBI Virus grows continuously, so a new download will not reproduce this
dataset exactly; the analysed file is therefore included.

---

## Software

Python 3.12 with `numpy`, `pandas`, `matplotlib` and `biopython`. No external
tools are required.

## Running

From inside this folder:

```bash
python a29l_haplotypes_web.py
```

The script reads `sequences.fasta`, writes all results to `results_web/`, and
prints the complete report, which is also saved as `results_web/qc_report.txt`.
It runs in a few seconds. Rerunning the script on the included data reproduces
every published number exactly.

---

## Method

**Reference.** Sequences are compared with the MPXV A29 reference sequence
YFT67316.1 (clade I; H74, R107). In the script, the reference sequence is
written out in full and labelled with the UniProt accession **Q77HM6**, which
has an identical sequence.

**Quality control.** Every count is reported before anything is removed, and
every excluded sequence is written to `excluded_sequences.csv` with its reason:

1. Sequences whose length differs from 110 aa are excluded.
2. Sequences containing non-standard or ambiguous residues are excluded.
3. Each remaining sequence is compared with the reference position by position
   over all 110 residues (ungapped identity; all sequences have the same
   length). Sequences below **90% identity** are excluded as unrelated proteins
   of coincidental length.

**Haplotypes and variable positions.** Haplotypes are defined by exact amino
acid sequence and numbered by frequency (H1 = most common). Differences are
reported relative to the dominant haplotype. A position is variable if more
than one residue occurs at it.

**Clade assignment.** Clade labels (I, Ia, Ib, II, IIa, IIb) are taken from the
NCBI Virus headers and grouped into clades I and II; the nomenclature
coincides with Nextstrain. Because residues 74 and 107 separate the two clades
without exception in labelled sequences, clade is imputed from this residue
pair for sequences lacking a label. Imputed clades are used **only for
denominators** in the temporal analysis and are always reported as imputed.

**Submission-batch diagnostic.** For each low-frequency haplotype with at least
five sequences (excluding the two clade-defining haplotypes), the script
checks whether it comes from a single GenBank deposit, using accession prefixes
and the span of accession numbers as a proxy for submission batches. A variant
confined to one deposit may reflect a single study rather than spread.

**Temporal frequency.** For the same haplotypes, the frequency per collection
year is computed, with the denominator being all sequences of the same
(imputed) clade in that year.

> Haplotype frequencies measure sequencing effort, not prevalence: they
> reflect which outbreaks and regions were sequenced and deposited.

Domain boundaries used for annotation: NTR 1–20, HBD 21–34, spacer 35–42, CCD
43–84, LZD 85–110.

---

## Results

| Step | Sequences |
| --- | --- |
| Downloaded (all 110 aa) | 9,870 |
| Excluded: ambiguous residue (B) | 1 |
| Excluded: < 90% identity to the reference | 184 |
| **Final curated dataset** | **9,685** |

Clade labels in the curated set: IIb 8,347; Ia 455; Ib 152; IIa 25; I 14;
II 9; unassigned 683.

**17 variable positions** define **17 haplotypes**:

| Haplotype | n | % | Difference from dominant | Clade |
| --- | --- | --- | --- | --- |
| H1 | 8,816 | 91.03 | dominant (R74, H107) | II |
| H2 | 716 | 7.39 | R74H; H107R | I |
| H3 | 104 | 1.07 | A91V | II |
| H4 | 19 | 0.20 | S21F | II |
| H5 | 9 | 0.09 | R74H; H107R; P108L | I |
| H6–H17 | 1–3 each | | single substitutions or pairs; 6 singletons | |

The complete table is `results_web/haplotypes.csv`, and every variable position
is listed in `results_web/variable_positions.csv`.

**Residues 74 and 107 discriminate clade perfectly.** Among the 9,002
sequences with a clade label, all 621 clade I sequences carry H74/R107 and all
8,381 clade II sequences carry R74/H107, with **zero discordant sequences**. The
683 unassigned sequences divide into 576 R74/H107 (imputed clade II) and 107
H74/R107 (imputed clade I), giving 8,957 clade II and 728 clade I sequences in
total.

**Low-frequency variants.** A91V (H3, 104 sequences) occurs across 24
separate deposits and three collection years (2024–2026), so it is not the
artefact of a single submission; S21F (H4, 19 sequences) spans six deposits,
all from 2024. H5 comes from a single deposit and should be treated with
caution. Details are in `batch_diagnostic.csv` and `temporal_frequencies.csv`.

---

## Output files (`results_web/`)

| File | Content |
| --- | --- |
| `qc_report.txt` | Complete printed report of the run |
| `excluded_sequences.csv` | The 185 excluded sequences, with reasons |
| `per_sequence.csv` | Every curated sequence: accession, clade, date, haplotype, residues 74/107, imputed clade |
| `haplotypes.csv` | All 17 haplotypes with counts and differences |
| `variable_positions.csv` | All 17 variable positions, with domains |
| `haplotype_by_clade.csv` | Haplotype × clade label counts |
| `batch_diagnostic.csv` | Submission-batch diagnostic |
| `temporal_frequencies.csv` | Frequency per year with clade denominators |
| `a29l_haplotype_analysis.png/.pdf/.svg` | Four-panel summary figure |

## Other files

`haplotypes17.fa` contains the sequences of the 17 haplotypes, prepared
manually in Jalview from the survey results, to display all haplotypes in a
single alignment.
