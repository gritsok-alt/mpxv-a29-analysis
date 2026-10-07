# Methods Summary — A29L Ortholog Discovery, Alignment and Effective Depth

A record of what each phase does, why, and what its limitations are. It
supports sections 2.1.1–2.1.2 of the article. Parameter values are those set
in each notebook's `Config`.

---

## Overview

The pipeline identifies orthologs of the monkeypox virus **A29L** protein
across *Poxviridae*, curates them, reduces redundancy, aligns them, and
characterises their conservation and effective alignment depth. A29L is the
MPXV ortholog of vaccinia **A27**, a short (~110 aa) IMV membrane protein
comprising an N-terminal heparin/GAG-binding region, a coiled-coil
oligomerisation domain, and a C-terminal leucine zipper that anchors it to the
virion through A17.

The central analytical challenge is discriminating true A29L orthologs from two
confounders: the **A26** paralog (a much larger protein that forms a
disulfide-linked complex with A27, and whose C-terminal region also carries a
PF02346 domain) and fragments of the unrelated **A-type inclusion (ATI)**
protein. Both are addressed with curated Pfam profiles.

Residue numbering throughout refers to the MPXV reference **YFT67316.1**.

---

## Phase 1 — Reference database construction

A non-redundant protein database is assembled for *Chordopoxvirinae* (NCBI taxid
`10241`, 266 species-level taxa) from two independent sources, both queried on
**22 July 2026**:

- **NCBI Protein**, enumerated per species. Species-level taxa are obtained
  from the taxonomy subtree, then each is searched separately for records
  within a **60–300 aa** window.
- **UniProtKB**, retrieved in one streamed query over the same taxon and
  length window.

Records are merged by **SHA-256 hash of the uppercased sequence**, so
identical sequences reported by both sources collapse into a single entry that
retains every contributing accession and organism name. Sequences containing
characters outside the standard amino-acid alphabet (plus the accepted
ambiguity codes `X B Z J U O`) are excluded. The resulting database contains
**10,017 sequences**.

**Per-species cap.** Retrieval is capped at **400 records per species**.
Heavily sequenced species (vaccinia, variola) otherwise contribute thousands of
near-identical strain deposits, many of them partial or misannotated. Greater
sampling depth yielded predominantly redundant or incorrectly annotated
sequences without identifying additional A29L orthologs.

**Length window.** The 60–300 aa window encompasses all known A29L homologues
(74–237 residues) while excluding the substantially longer A26L paralogues.

> **Limitation.** Fragments of a ~110 aa protein pass easily through the
> length window; it is the per-species cap that removes most fragmentary
> records. Selection within a species follows NCBI's own result ordering.

**Outputs:** `search_db.fasta` (non-redundant database), `accession_map.csv`
(provenance: hash, all contributing accessions, organisms), and the raw source
snapshots `genbank_perspecies.fasta` and `uniprot_chordo.fasta`.

---

## Phase 2 — Homology search

The A29L query is searched against the Phase 1 database by three methods with
complementary sensitivity profiles:

| Method | Version | Settings | Rationale |
| --- | --- | --- | --- |
| PSI-BLAST | BLAST+ 2.16.0 | 5 iterations; E and inclusion ≤ 1e-3; BLOSUM62; gap open 11, extend 1 | Detects divergent homologs via a query-derived profile |
| jackhmmer | HMMER 3.4 | 5 iterations; E and incE ≤ 1e-3 | Different sensitivity; catches hits PSI-BLAST misses |
| hmmsearch vs PF02346 | HMMER 3.4 | E ≤ 1e-3 | Curated profile, independent of the query sequence |

Hits are unioned into one candidate set, and `search_hits_provenance.csv`
records which method(s) found each sequence. All three methods converged on the
same **82 candidate sequences**.

---

## Phase 3 — Search of unassigned taxa

*Chordopoxvirinae* and *Entomopoxvirinae* do not exhaust *Poxviridae* (taxid
`10240`); some species are assigned to neither. These taxa are computed as

> Poxviridae species − Chordopoxvirinae species − Entomopoxvirinae species

giving **40 unassigned taxa**, whose proteins are retrieved and searched with
jackhmmer under the same thresholds. This ensures that a genuine ortholog in an
unassigned lineage is not missed merely because it falls outside the two named
subfamilies. The search identified **3 additional candidates**
(`orphan_hits.fasta`).

---

## Phase 4 — Unified domain classification

Candidates from Phases 2 and 3 are **pooled** and classified together, so hits
from unassigned taxa receive the same screening as the main set. Each
candidate is scanned with `hmmsearch --cut_ga` against three Pfam profiles
(downloaded from InterPro on 12 February 2026), and filters are applied in
order:

1. **PF06086** present (A26L/A30L family) → A26 paralog.
2. **PF04508** present (*Pox_A_type_inc*, ATI repeat) → ATI protein.
3. **Length > 250 aa** → A26 paralog. A29L is ~110 aa; VACV A26 is ~58 kDa
   versus A27's ~12 kDa.
4. Otherwise, carrying **PF02346** (*Vac_Fusion*) → **A29L ortholog**.

Orthology therefore requires PF02346 together with the absence of PF06086 and
PF04508. PF02346 alone is insufficient, because the C-terminal region of A26L
also contains this domain.

Of the **85 pooled candidates** (82 from Phase 2 and 3 from Phase 3), **2 were
classified as A26 paralogs**, leaving **83 orthologs** (`true_orthologs.fasta`).

The A29L C-terminal motif (`[KR]?KID[VM]QTG`) is recorded as a per-sequence
annotation but **is not used for classification**.

> **Note.** The motif was derived in-house from the query alignment and
> corresponds to no published or registered signature. Keeping it as an
> annotation means no classification decision rests on an unpublished
> pattern. The Pfam accessions, by contrast, are citeable and stable.

> **Limitation.** The Phase 1 length window already excludes full-length A26
> (~500 aa) and ATI (~700+ aa), so these filters act as a safety net catching a
> small number of residual fragments rather than performing the bulk of the
> discrimination. `candidates_classified.json` records the class assigned to
> every candidate.

---

## Phase 5 — Curation

Curation only; nothing is clustered or reordered. Three steps:

**Taxonomic annotation.** Each accession's organism is resolved in tiers (the
Phase 1 accession map, then the NCBI record, then UniProt) and its genus taken
from the NCBI taxonomy lineage. Headers are rewritten as
`Genus|Accession|Species`. Taxa with **no genus rank in NCBI** keep their
species and are labelled `unclassified` rather than being forced to a higher
rank. Assignments are then **hand-verified** through an editable table
(`taxonomy_overrides.csv`) that records both NCBI's value and the curated one,
and naming is standardised to ICTV binomial style.

**Domain QC.** Every sequence is scored against PF02346 and PF06086, each
scanned twice: with `--cut_ga` (Pfam gathering threshold) for pass/fail, and
without a cutoff so that below-threshold sequences still report a bit score.
`qc_scores.csv` records both. **Nothing is removed at this step**; it flags
candidates for review.

**Manual exclusion.** All removals happen here, through `manual_exclusions.py`,
each with a documented reason. Of the 83 orthologs, **nine sequences were removed**: four with
anomalous lengths, one containing a C-terminal frameshift, one with ten
ambiguous residues, one lacking part of the N-terminus, and two
*Crocodylidpoxvirus* proteins below the PF02346 gathering threshold, which
therefore could not be confidently assigned to the family. The final curated
set contains **74 sequences from 43 species**.

> A true A29L ortholog is expected to pass PF02346 and fail PF06086. The
> PF06086 scan is independent evidence against A26 paralogy; the PF02346 scan
> *corroborates* rather than independently confirms, since the same family
> defined the ortholog set in Phase 4.

---

## Phase 6 — Redundancy reduction

The curated set is clustered with **CD-HIT 4.8.1** at 99, 98, 95 and 90%
identity (`-n 5 -d 0 -M 2000`; word size 5 is valid for all protein thresholds
≥ 0.7). Producing all four supports a redundancy-sensitivity check: a
downstream result stable across thresholds is not an artefact of sequence
redundancy.

**Reference hold-out.** The two MPXV reference sequences (`YFT67316.1`,
`YFT66960.1`) are removed before clustering and reinserted into every
non-redundant set, guaranteeing their presence in all downstream analyses.

**Species-coverage verification.** CD-HIT clusters by identity, not taxonomy,
so distinct species with near-identical A29L can collapse into one
representative. Coverage is therefore reported for each threshold.

Result (curated input: 74 sequences spanning 43 species):

| Threshold | Representatives | Reduction | Species covered |
| --- | --- | --- | --- |
| 99% | 65 | 12.2% | 43/43 |
| **98%** | **60** | **18.9%** | **43/43** |
| 95% | 45 | 39.2% | 34/43 |
| 90% | 38 | 48.6% | 32/43 |

Representatives include the two held-out references; reduction is calculated
as 1 − representatives / 74.

> **Note.** The reduction percentages printed in `output_phase6/phase6.log`
> (14.9, 21.6, 41.9 and 51.4%) divide the number of clusters *excluding* the
> held-out references by the total *including* them. The values above use the
> consistent calculation. No result depends on this summary line.

**98% was selected** as the most stringent threshold at which no species is
lost. At 95%, nine species lose their own representative: *Orthopoxvirus*
abatino, akhmetapox, cowpox, ectromelia, taterapox, vaccinia and variola, and
the unclassified moosepox and white-tailed deer poxviruses. At 90%, goatpox and
lumpy skin disease viruses (*Capripoxvirus*) are lost in addition. This is a
direct consequence of A29L's high conservation within *Orthopoxvirus*.

---

## Phase 7 — Alignment and conservation analysis

The 98% representative set (60 sequences) is aligned with **MAFFT 7.525
L-INS-i** (`--localpair --maxiterate 1000`), the accuracy-oriented mode
appropriate at this scale. The resulting **300-column alignment is left
untrimmed**.

**Pairwise identity** is computed over columns where both sequences have a
residue, so gap positions do not dilute the score (`pairwise_identity.csv`).

**Genus heatmap.** The identity matrix is ordered by genus with block
boundaries drawn, separating within-genus conservation from between-genus
divergence.

**Conservation profiles.** For each column, the conservation score is the
frequency of the most common residue over the total number of sequences
(identity × occupancy), so a column scores high only if it is both
well-populated and uniform. Occupancy is overlaid. The profile is mapped onto
the MPXV reference, located **by sequence match rather than accession**, and
annotated with the historically annotated signal peptide, GAG-binding region,
spacer, coiled-coil domain, leucine zipper, and the disulfide-linked motif
(`A29L_reference_mapped_conservation.csv`).

Figures are exported as PNG at 1000 dpi on white backgrounds. The article
figures were redrawn from these outputs.

**Domain-level summary.** The per-column profiles are summarised over the
annotated structural domains to give Tables S4 and S5. For each domain the mean
occupancy, conservation score and identity where present are taken across its
residue positions, referenced to MPXV A29 (`YFT67316.1`). The summary is computed
at two levels: the full 60-sequence *Chordopoxvirinae* alignment (Table S4) and a
14-sequence *Orthopoxvirus* subset (Table S5), the latter comprising every
sequence assigned to genus *Orthopoxvirus* in the 98% representative set. The
subset is extracted as a **column subset of the 60-sequence alignment, without
realignment**, so that column numbering — and therefore residue correspondence —
is identical at both levels.

The calculation lives in `occupancy_conservation_identity/`, which holds
`occupancy_conservation_identity.py` together with its two inputs,
`60sequences_subfamily.fasta` and `14sequences_genus.fasta`. The script depends
only on the Python standard library. It refuses any input whose records differ in
length, and warns if the ungapped reference is not 110 residues — the check that
the subset was extracted correctly and that the reference survived into it.

> **Limitation.** The conservation score is identity weighted by occupancy, so a
> low value may reflect poor column coverage rather than sequence divergence. The
> two are separated by reading it against identity where present: at the
> *Chordopoxvirinae* level the GAG-binding domain falls from 0.457 to 0.290 through
> gaps alone. Fig. 2 plots the weighted score and should be read accordingly.

---

## Effective alignment depth (`neff/`)

The effective number of sequences is

> N_eff = Σ_i 1 / |{ j : I(i, j) ≥ 0.62 }|

where I(i, j) is the pairwise identity over mutually ungapped positions (X
counted as a gap), and each sequence counts as its own neighbour. The 62%
threshold follows PSICOV. N_eff/L uses L = 110, the length of the reference.

**A single alignment.** All values are computed from the Phase 7 alignment,
so that the only variable between rows is which sequences are present, never
where the gaps lie:

- **nr90, nr95, nr98** are strict subsets of the alignment, obtained by
  deleting rows.
- **nr99** is the alignment minus one sequence (`P20535.1`, vaccinia, for which
  CD-HIT chose a different representative at 99%), plus six sequences threaded
  in with `mafft --add`, giving 65 sequences.
- **nr100** is the full curated set: the 60 aligned sequences plus the 14
  remaining ones threaded in with `mafft --add`, giving 74 sequences.

```
mafft --add <extra.fasta> --localpair --maxiterate 1000 --amino --quiet <base.fasta>
```

Columns inserted by `--add` are gaps in all original sequences and are
therefore skipped when identity is computed between any two of them; the
published values for nr90, nr95 and nr98 are unchanged by construction.

**Core.** N_eff is also computed for the confidently aligned core, the columns
holding reference residues 44–110.

Published values: N_eff = 18.01 (nr90, 38 sequences), 17.91 (nr95, 45),
17.85 (nr98, 60) and 13.89 (nr98 core 44–110). All rows, including nr99 and
nr100, are listed in `neff/neff_table.csv`.

---

## Summary of documented limitations

1. Retrieval is capped at 400 records per species; this, not the length
   window, is what removes most fragmentary entries.
2. The A29L motif regex is an in-house pattern and is used only as an
   annotation.
3. The 250 aa length threshold is reasoned but not externally derived.
4. The PF02346 QC corroborates rather than independently validates ortholog
   assignment.
5. Genus assignment depends on NCBI taxonomy, which is mid-migration to ICTV
   binomials; assignments were hand-verified and standardised, with the raw
   NCBI values retained.
6. Species with no genus rank in NCBI are labelled `unclassified` rather than
   assigned a higher rank.
7. Clustering at ≤ 95% identity collapses distinct orthopoxvirus species; 98%
   was selected on the basis of the coverage check.
8. The subset and `base`/`toadd` files in `neff/` were prepared from the
   Phase 6 and Phase 7 outputs outside the notebooks; their composition is
   described above.
9. The conservation score plotted in Fig. 2 and reported in Tables S4–S5 is
   identity weighted by occupancy; low regions can arise from coverage rather than
   divergence, and identity where present is reported alongside to separate them.

---

## Reproducibility

Every phase is checkpointed and driven by a single `Config` object, so
parameters can be quoted directly from the notebooks. All curation decisions
(exclusions, renames, taxonomy corrections) are held in editable files
(`manual_exclusions.py`, `output_phase5/taxonomy_overrides.csv`) rather than in
code, providing an auditable record of every human judgement applied to the
dataset.

| Item | Value |
| --- | --- |
| Operating environment | Windows 11; Ubuntu 26.04 under WSL2 |
| Python | 3.12 |
| BLAST+ | 2.16.0 |
| HMMER | 3.4 (Aug 2023) |
| CD-HIT | 4.8.1 |
| MAFFT | 7.525 |
| Pfam profiles | PF02346, PF06086, PF04508; downloaded from InterPro, 12 Feb 2026 |
| NCBI / UniProt retrieval | 22 July 2026 |
