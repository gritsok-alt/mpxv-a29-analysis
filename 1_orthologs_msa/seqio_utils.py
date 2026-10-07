"""Shared sequence and search-output parsers.

Small, dependency-free helpers used by more than one phase, kept here
so every phase parses FASTA and HMMER/BLAST tabular output the same way.
"""

from __future__ import annotations

from pathlib import Path


def parse_fasta_records(path: Path) -> list[tuple[str, str]]:
    """Return each FASTA record as a (header, sequence) tuple."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_lines: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header, seq_lines = line[1:].strip(), []
            elif line.strip():
                seq_lines.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records


def parse_fasta_by_accession(path: Path) -> dict[str, str]:
    """Return an accession-keyed dict for O(1) lookup by first header token."""
    seqs: dict[str, str] = {}
    header: str | None = None
    for hdr, seq in parse_fasta_records(path):
        header = hdr.split()[0]
        seqs[header] = seq
    return seqs


def parse_hmmer_tabular(path: Path) -> set[str]:
    """Return target names from an HMMER --tblout file (jackhmmer/hmmsearch)."""
    hits: set[str] = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            hits.add(line.split()[0])
    return hits


def parse_domtblout(path: Path) -> set[str]:
    """Return target names from an HMMER --domtblout file.

    The target name is the first whitespace-delimited field; '#'-comment
    and blank lines are skipped.
    """
    hits: set[str] = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            hits.add(line.split()[0])
    return hits


def parse_blast_tabular(path: Path) -> set[str]:
    """Return unique subject accessions from a BLAST -outfmt 6 file.

    Tolerates PSI-BLAST's per-iteration blank lines, '#'-comment lines,
    and the 'Search has CONVERGED!' marker written between rounds.
    """
    hits: set[str] = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line == "Search has CONVERGED!":
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                hits.add(parts[1])
    return hits


def count_cdhit_clusters(clstr_path: Path) -> int:
    """Return the number of clusters in a CD-HIT .clstr file.

    Each cluster begins with a '>Cluster' header line; counting those gives
    the number of representative sequences produced.
    """
    with open(clstr_path, encoding="utf-8") as fh:
        return sum(1 for line in fh if line.startswith(">Cluster"))


def load_accession_organisms(accession_map_csv: Path) -> dict[str, set[str]]:
    """Map every accession in a Phase 1 accession_map.csv to its organism set.

    Each row lists a preferred accession plus semicolon-separated ncbi_accs and
    uniprot_accs, all referring to the same merged sequence, and a
    ' | '-separated organisms field. Every one of those accessions is indexed
    to that row's organism set, so a lookup succeeds whichever accession form a
    downstream file carries.
    """
    import csv

    acc_to_orgs: dict[str, set[str]] = {}
    with open(accession_map_csv, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            orgs = {o.strip() for o in row.get("organisms", "").split("|") if o.strip()}
            accs = {row.get("preferred_acc", "").strip()}
            for field_name in ("ncbi_accs", "uniprot_accs"):
                accs.update(a.strip() for a in row.get(field_name, "").split(";") if a.strip())
            for acc in accs:
                if acc:
                    acc_to_orgs.setdefault(acc, set()).update(orgs)
    return acc_to_orgs


def parse_domtbl_scores(path: Path) -> dict[str, float]:
    """Map each target to its best full-sequence bit score from a --domtblout file.

    In hmmsearch --domtblout output the target name is column 1 and the
    full-sequence score is column 8 (1-based). A target may appear on several
    domain rows; the highest full-sequence score seen is kept. '#'-comment and
    blank lines are skipped.
    """
    best: dict[str, float] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 8:
                continue
            target = cols[0]
            try:
                score = float(cols[7])
            except ValueError:
                continue
            if target not in best or score > best[target]:
                best[target] = score
    return best


def clean_accession(raw: str) -> str:
    """Reduce a UniProt/GenBank-style pipe header to its bare accession.

    'sp|P26654.1|VFUS_ORFN2' -> 'P26654.1'; 'tr|X|Y' and 'db|X|Y' likewise.
    A token with no recognised database prefix is returned unchanged, so bare
    accessions pass through untouched.
    """
    token = raw.split()[0]
    if "|" in token:
        parts = token.split("|")
        if len(parts) >= 3 and parts[0].lower() in {"sp", "tr", "db", "gb", "ref", "pdb", "emb", "dbj"}:
            return parts[1]
    return token
