# Step 1 — A29L ortholog discovery, curation and alignment

This folder contains the code, curation decisions and results for sections
2.1.1–2.1.2 of the article: construction of a poxvirus protein database,
identification and curation of orthologs of the monkeypox virus **A29L**
protein (vaccinia A27 family), redundancy reduction, multiple sequence
alignment, conservation analysis, and effective alignment depth (N_eff).

The central output of this step, used by all later steps, is the final
alignment **`output_phase7/A29L_msa_linsi.fasta`** (60 sequences × 300
columns, MAFFT L-INS-i, untrimmed).

A detailed description of each phase, its rationale and its limitations is
given in [`METHODS_SUMMARY.md`](METHODS_SUMMARY.md).

---

## Folder layout

```
1_orthologs_msa/
├── README.md                            this file
├── METHODS_SUMMARY.md                   methods, rationale and limitations
├── phase1_database.ipynb                build the protein database
├── phase2_homology_search.ipynb         PSI-BLAST, jackhmmer, hmmsearch
├── phase3_orphan_taxa.ipynb             search of unassigned Poxviridae taxa
├── phase4_domain_classification.ipynb   Pfam-based orthology assignment
├── phase5_curation.ipynb                taxonomy, domain QC, manual exclusion
├── phase6_clustering.ipynb              CD-HIT redundancy reduction
├── phase7_alignment_conservation.ipynb  MAFFT alignment and conservation
├── ncbi_fetch_utils.py                  shared NCBI/UniProt fetching
├── seqio_utils.py                       shared FASTA / HMMER / BLAST parsers
├── manual_exclusions.py                 all curation decisions (exclusions, renames)
├── PF02346.hmm                          Pfam: Vac_Fusion (A27/A29L family)
├── PF06086.hmm                          Pfam: A26L/A30L family (A26 paralog marker)
├── PF04508.hmm                          Pfam: Pox_A_type_inc (ATI repeat)
├── output_phase1/ … output_phase7/      results of each phase
├── neff/                                effective alignment depth (N_eff)
└── occupancy_conservation_identity/
    ├── occupancy_conservation_identity.py   domain-level summary (Tables S4, S5)
    ├── 60sequences_subfamily.fasta          Chordopoxvirinae alignment (60 seq)
    └── 14sequences_genus.fasta              Orthopoxvirus subset (14 seq)
```

---

## Software

The analysis was run on Windows 11 with the external tools installed in
Ubuntu 26.04 under WSL2.

| Software | Version | Used in |
| --- | --- | --- |
| Python | 3.12 | all notebooks and scripts |
| BLAST+ (`makeblastdb`, `psiblast`) | 2.16.0 | Phase 2 |
| HMMER (`jackhmmer`, `hmmsearch`) | 3.4 (Aug 2023) | Phases 2–5 |
| CD-HIT | 4.8.1 | Phase 6 |
| MAFFT | 7.525 | Phase 7, `neff/` |

Python packages:

```
pip install numpy pandas matplotlib nbformat biopython
```

The three `.py` modules must stay **in this folder, beside the notebooks**.
Open this folder as the working directory before running any notebook.

External tools are called through WSL via a `wsl_prefix` field in each
notebook's `Config`. On a native Linux system, set `wsl_prefix` to `[]`.
Missing tools can be installed inside WSL with
`sudo apt install ncbi-blast+ hmmer cd-hit mafft`.

### Pfam profiles

The three Pfam profiles are **included** in this folder. They were downloaded
from [InterPro](https://www.ebi.ac.uk/interpro/) on 12 February 2026. Pfam
accessions are stable, but family names can change between releases.

### NCBI credentials

Phases 1, 3 and 5 query NCBI. An NCBI API key
([account settings](https://www.ncbi.nlm.nih.gov/account/settings/)) raises
the rate limit from 3 to 10 requests per second. Credentials are read from
environment variables and are never stored in the code:

```bash
export NCBI_EMAIL="you@example.org"
export NCBI_API_KEY="your_key_here"
```

Alternatively, call `ncbi.configure(email=..., api_key=...)` in the notebook
before any fetch. `ncbi.rate_limit_s()` returns `0.11` when the key is found
and `0.34` when it is not. The pipeline also runs without a key, only slower.

---

## Data provenance

NCBI and UniProt were queried on **22 July 2026**. Because both databases
change continuously, a new run will not reproduce the retrieval exactly. The
raw snapshots of that retrieval are therefore included
(`output_phase1/genbank_perspecies.fasta`, `output_phase1/uniprot_chordo.fasta`,
`output_phase3/orphan_proteins.fasta`), together with the exact identifiers
and taxa queried (`all_ncbi_ids.json`, `sp_taxids.json`, `orphan_ids.json`,
`orphan_taxids.json`). All downstream phases can be rerun from these files.

The BLAST database files built by `makeblastdb` in `output_phase2`
(`searchdb.*`) are not included; they are regenerated when Phase 2 is run.

---

## Running the pipeline

Run the notebooks in order, each from top to bottom. Every phase reads its
inputs from the previous phase's `output_phaseN/` folder by relative path, so
no files need to be moved. Phase 5 includes a manual curation step (see below).

| Phase | Notebook                        | Key output                                                                                        |
|-------|---------------------------------|---------------------------------------------------------------------------------------------------|
| 1     | `phase1_database`               | `search_db.fasta` (10,017 sequences), `accession_map.csv`                                         |
| 2     | `phase2_homology_search`        | `homology_candidates.fasta`, `search_hits_provenance.csv`                                         |
| 3     | `phase3_orphan_taxa`            | `orphan_hits.fasta` (3 sequences)                                                                 |
| 4     | `phase4_domain_classification`  | `true_orthologs.fasta` (83 of 85 pooled candidates), `candidates_classified.json`                 |
| 5     | `phase5_curation`               | `manually_curated_orthologs.fasta` (74 sequences), `taxonomy_overrides.csv`, `qc_scores.csv`      |
| 6     | `phase6_clustering`             | `orthologs_nr99/98/95/90.fasta` and `.clstr` files                                                |
| 7     | `phase7_alignment_conservation` | `A29L_msa_linsi.fasta`, `pairwise_identity.csv`, `A29L_reference_mapped_conservation.csv`, figures |
| -     | -                               | Tables S4 and S5, optional per-position CSVs files|
The figures produced by Phase 7 (`A29L_identity_heatmap.png`,
`A29L_global_conservation_profile.png`, `A29L_domain_conservation_profile.png`)
are the basis of the corresponding article figures, which were redrawn for
publication. The underlying values are in the accompanying CSV files.

### Manual curation in Phase 5

Phase 5 is a two-pass, human-in-the-loop phase:

1. The **first run** writes `output_phase5/taxonomy_overrides.csv` (a
   pre-filled taxonomy table) and `qc_scores.csv` (per-sequence Pfam bit
   scores).
2. **Review both.** In `taxonomy_overrides.csv`, edit only the `genus`,
   `species` and `notes` columns. Never edit `accession` (the join key), and
   leave `ncbi_genus`, `ncbi_species` and `provenance` as the record of what
   NCBI reported. Save as UTF-8 CSV. In `qc_scores.csv`, the `flag` values
   `PF02346_below_GA` and `PF06086_hit` mark sequences worth inspecting.
3. **Add sequences to be removed** to `manual_exclusions.py`, each with a
   documented reason.
4. **Rerun.** Edits are applied and never overwritten.

Delete `taxonomy_overrides.csv` to rebuild the pre-filled table from scratch.

### Editing `manual_exclusions.py`

This file holds every curation decision and is the only place where
sequences are removed or relabelled.

```python
EXCLUSIONS = {
    "QGT49194.1": "Crocodylidpoxvirus, below PF02346 GA threshold",
}

RENAMES = {}
```

Every entry in `EXCLUSIONS` must carry a reason. `RENAMES` relabels an
accession; use it sparingly, since it changes the stated provenance of a
sequence, and only after verifying that both records hold an identical
sequence. Malformed UniProt headers (`sp|P26654.1|VFUS_ORFN2`) are cleaned to
their bare accession automatically.

### Domain-level summary (Tables S4 and S5)

Run from inside `occupancy_conservation_identity/`:

```bash
python3 occupancy_conservation_identity.py
```

Both alignments are found by their default filenames and both tables are printed.
Add `--per-position` for the underlying vectors, `--csv` to write four CSVs beside
the script, or `--ref ACCESSION` if the reference header changes. Two file paths
given as arguments override the defaults.

No external tools are required and NumPy is not needed — the script uses only the
Python standard library, so it runs without the rest of the pipeline's
environment.

---

## Effective alignment depth (`neff/`)

`neff/compute_neff.py` computes N_eff for the alignment at each CD-HIT
threshold, and for the confidently aligned core (reference residues 44–110).
See `neff/README.md` for details. All values are computed from the single
Phase 7 alignment: smaller sets by deleting rows, larger sets by threading the
additional sequences into it with `mafft --add`, so that alignment columns
never change between rows.

To reproduce the table, run from inside `neff/`:

```bash
python compute_neff.py
```

It prints the table and writes `neff_table.csv`. The published values are
N_eff = 18.01 (nr90), 17.91 (nr95), 17.85 (nr98) and 13.89 (nr98, core 44–110).

---

## Configuration

| Setting | Phase | Value | Meaning |
| --- | --- | --- | --- |
| `max_ids_per_species` | 1, 3 | `400` | Cap on records fetched per species |
| `min_len` / `max_len` | 1, 3 | `60` / `300` | Length window at retrieval |
| `length_threshold` | 4 | `250` | Above this, classified as A26 paralog |
| `hold_out_accessions` | 6 | `YFT67316.1`, `YFT66960.1` | MPXV references, never clustered away |
| `cdhit_identities` | 6 | `0.99, 0.98, 0.95, 0.90` | Clustering thresholds |
| `reference_seq` | 7 | MPXV A29L (`YFT67316.1`) | Reference, matched by sequence |
| `figure_dpi` | 7 | `1000` | PNG resolution |

---

## Checkpointing and reruns

Each phase skips steps whose output already exists, so an interrupted run
resumes rather than restarting. To force a step to rerun, delete its output:

| To redo | Delete |
| --- | --- |
| NCBI ID discovery | `output_phase1/all_ncbi_ids.json` |
| Taxonomy pre-fill | `output_phase5/taxonomy_overrides.csv` |
| Clustering | `output_phase6/` |
| Alignment | `output_phase7/A29L_msa_linsi.fasta` |

Deleting a phase's whole output folder forces a clean rerun of that phase.

---

## Troubleshooting

**Figures have no visible text.** A dark matplotlib style applied globally
sets text to white, which is invisible on the white publication background.
Phase 7's style cell calls `matplotlib.rcdefaults()` and forces black text;
run the notebook from the top so that this cell executes before the figure
cells. `print(plt.rcParams["text.color"])` must print `black`.

**`Reference A29L sequence not found in the alignment`.** The reference is
matched by sequence, not accession. Confirm that it is present in
`orthologs_nr98.fasta`; Phase 6 logs the line "Holding 2 reference
sequence(s) out of clustering".

**Genus shows `unclassified`.** Either NCBI has no genus rank for the taxon,
or the organism name could not be resolved; the `provenance` column
distinguishes the two. Correct it by hand in `taxonomy_overrides.csv`.

**Species lost during clustering.** Expected below 98% identity, where
closely related orthopoxviruses collapse into shared representatives. Phase 6
lists the affected species for each threshold.

**A tool "is not recognised".** The tool is missing inside WSL, not Windows.
Install it in the WSL distribution with `apt`.
