# N_eff — one number per CD-HIT threshold

Goal: six rows in a table. Every row must come from **the same alignment**, so
that the only thing changing between rows is *which sequences are present* —
never where the gaps sit.

The alignment is `MSA/output_phase7/A29L_msa_linsi.fasta` (60 sequences ×
300 columns, MAFFT L-INS-i, the nr98 set). Nothing in this folder re-aligns it.

---

## What is already done

Three of the five sets are just **rows taken out of** that alignment, because
nr90 and nr95 are strict subsets of nr98. No software needed — already written:

| file | sequences |
|---|---|
| `aln_nr90.fasta` | 38 |
| `aln_nr95.fasta` | 45 |
| `aln_nr98.fasta` | 60 (the original, unchanged) |

## What you have to run — two commands

nr99 and nr100 contain sequences that are **not** in the alignment, so they
cannot be made by deleting rows. They are made by *threading the extra sequences
into the existing alignment* with `mafft --add`, which leaves all 300 original
columns exactly as they are.

Open WSL in this folder and run:

    mafft --add nr100_toadd.fasta --localpair --maxiterate 1000 --amino --quiet nr100_base.fasta > aln_nr100.fasta

    mafft --add nr99_toadd.fasta  --localpair --maxiterate 1000 --amino --quiet nr99_base.fasta  > aln_nr99.fasta

Then check the counts are 74 and 65:

    grep -c ">" aln_nr100.fasta aln_nr99.fasta

If MAFFT refuses `--localpair` together with `--add`, delete
`--localpair --maxiterate 1000` and run it again. It changes nothing that
matters here.

Do **not** add `--keeplength`. It would force the output back to 300 columns by
deleting residues from the sequences you just added.

## Then

    python compute_neff.py

It prints the table and writes `neff_table.csv`. It skips any file that isn't
there, so you can run it now, before the MAFFT step, to see the first four rows.

---

## Why the two `_base` files are different sizes

`nr100_base.fasta` is the alignment untouched (60). All 60 are in the 74-sequence
curated set, so nr100 = alignment + 14 added.

`nr99_base.fasta` is the alignment **minus one sequence** (59). CD-HIT at 99 %
picked a different representative for one vaccinia cluster, so
`Orthopoxvirus|P20535.1|Orthopoxvirus_vaccinia` is in nr98 but not in nr99.
It has already been removed for you. nr99 = 59 + 6 added = 65.

## Why not just re-align each set separately?

A smaller set would probably align slightly *better*. That is the problem, not
the benefit. The claim being tested is "N_eff does not depend on the clustering
threshold". Re-aligning changes two things at once — which sequences are present
*and* where the gaps go — so any flatness or slope in the table becomes
uninterpretable. Holding the alignment fixed leaves membership as the only
variable, which is the same reason you don't re-run an analysis at each step of
a rarefaction curve.

## Why adding sequences cannot change the numbers you already published

Identity is measured only over columns where **both** sequences have a residue.
`--add` may insert new columns, but those columns are gaps in all 60 original
sequences, so every pair of original sequences skips them. Their pairwise
identities are therefore unchanged, and so is N_eff for nr90, nr95 and nr98.
`compute_neff.py` prints the expected values so you can confirm this.
