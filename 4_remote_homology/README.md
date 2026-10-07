# Step 4 — Remote homology search (HHblits and HHpred)

This folder contains the results of the remote-homology search described in
section 2.1.5 of the article, which asked whether A29L has detectable
homologues outside the *Poxviridae*.

Both searches were run on the **MPI Bioinformatics Toolkit** web service on
**28 July 2026**, with the MPXV A29L reference **YFT67316.1** (110 aa) as query.
There is therefore no code in this folder. Toolkit jobs are deleted after some
time, so the downloaded result files below are the permanent record of the
analysis. The job numbers identify the original Toolkit runs.

---

## Folder layout

```
4_remote_homology/
├── README.md
├── hhblits/                          job 9121590
│   ├── hhblits_9121590.out           full search report
│   ├── hhblits_full_9121590.a3m      collected homologues (A3M alignment)
│   └── hhblits_fullQT_9121590.a3m    query-template version of the alignment
└── hhpred/                           job 8650667
    ├── hhpred_8650667.hhr            template hits with probabilities and E-values
    ├── hhpred_full_8650667.a3m       query alignment used for the search
    └── hhpred_querytemplate_8650667.fasta   query–template alignments
```

The `.hhr`, `.out` and `.fasta` files are plain text; `.a3m` alignments can be
opened in Jalview or any text editor.

---

## Searches and settings

**1. HHblits against UniRef30** (job 9121590) — iterative profile search for
homologous sequences.

| Setting | Value |
| --- | --- |
| Database | UniRef30_2023_02 |
| Iterations | 3 (`-n 3`) |
| E-value inclusion threshold | 1e-3 (`-e 1e-3`) |
| Other parameters | `-p 10 -Z 1000 -z 1 -b 1 -B 1000` |

**2. HHpred / HHsearch against structure and domain databases** (job 8650667).
On the Toolkit, HHpred builds the query alignment with HHblits iterations and
then runs a single HHsearch profile–profile comparison against the template
databases.

| Setting | Value |
| --- | --- |
| Databases | PDB_mmCIF70_20_Feb, Pfam-A v38.2, NCBI Conserved Domains v3.19 |
| Parameters | `-p 10 -Z 1000 -loc -z 1 -b 1 -B 1000 -ssm 0 -sc 1 -seq 1 -norealign -maxres 32000` |
| Query alignment | 60 of 62 sequences, Neff = 3.26 |
| HMMs searched | 123,631 |

The complete HHsearch command line is recorded at the top of
`hhpred_8650667.hhr`.

---

## Results

**HHblits.** Three iterations against UniRef30 produced a profile of only 55
sequences (55 of 63 retained; Neff = 3.21). Of the 1,000 hits reported,
26 were significant (E ≤ 1e-3), and **all 26 are poxvirus proteins**: IMV
envelope/fusion proteins (A27/A29L orthologs such as vaccinia OPG155,
molluscum contagiosum MC131/MC133) and A-type inclusion (P4c) proteins, from
orthopoxviruses, parapoxviruses, capripoxviruses, avipoxviruses and other
chordopoxviruses. The best-scoring hits outside the *Poxviridae* (uncharacterised
proteins, a haemolysin, DUF1664 and SlyX-family proteins) have E-values of
0.3 or higher.

**HHpred.** HHsearch returned 526 hits, of which only two were significant
(E ≤ 1e-3), both expected:

| Rank | Hit | Probability (%) | E-value | Query residues | Template residues |
| --- | --- | --- | --- | --- | --- |
| 1 | PF02346 (Vac_Fusion, Pfam) | 99.7 | 6.5e-21 | 50–105 | 1–49 of 49 |
| 2 | 3VOP chain C (VACV A27, PDB) | 99.4 | 6.2e-17 | 21–84 | 1–64 of 64 |
| 3 | 3EFG chain A (SlyX homolog, PDB) | 88.1 | 0.38 | 44–97 | 9–62 of 78 |

The remaining 524 hits were predominantly generic coiled-coil proteins (SlyX,
ZapB, TRAF coiled-coil domains, designed coiled-coil trimers, nucleoporins),
frequently aligned to query residues ~45–100. Their alignment reflects shared
helical, coiled-coil architecture rather than homology.

**Conclusion.** No homology to A29L was detectable outside the *Poxviridae*.
At the sensitivity of current profile–profile searches, A29L appears to be an
orphan protein, and the only homologous structural template is the vaccinia
A27 structure (3VOP).

> **Note on significance.** Hits are called significant here by E-value
> (≤ 1e-3). HHpred users often also consult the probability column, in which
> the best coiled-coil hit (3EFG, SlyX homolog) reaches 88.1%, at an E-value of
> 0.38. High probabilities for short coiled-coil matches are a known feature of
> HHpred and reflect structural similarity of heptad repeats rather than
> common ancestry.
