#!/usr/bin/env python3
"""Tables S4 and S5 — occupancy and conservation by structural domain.

Computes, for each alignment, the per-position occupancy, conservation score and
identity where present, referenced to MPXV A29 (YFT67316.1), and reports the mean
of each over the structural domains.

Default inputs, expected in the same directory as this script:
    60sequences_subfamily.fasta   -> Table S4 (Chordopoxvirinae level)
    14sequences_genus.fasta       -> Table S5 (Orthopoxvirus level)

Usage:
    python3 occupancy_conservation_identity.py
    python3 occupancy_conservation_identity.py --per-position
    python3 occupancy_conservation_identity.py --csv
    python3 occupancy_conservation_identity.py --ref OTHER_ACCESSION
    python3 occupancy_conservation_identity.py aln1.fasta aln2.fasta

Prints to stdout. Writes files only with --csv. Standard library only.
"""

import sys
from collections import Counter
from pathlib import Path

GAPS = set("-.~")
DEFAULT_REF = "YFT67316.1"

DATASETS = [
    ("Table S4", "60sequences_subfamily.fasta", "Chordopoxvirinae, 60 orthologues"),
    ("Table S5", "14sequences_genus.fasta", "Orthopoxvirus, 14 orthologues"),
]

DOMAINS = [
    ("N-terminal region",      1,  20),
    ("GAG-binding domain",     21, 32),
    ("Spacer",                 33, 42),
    ("Coiled-coil domain",     43, 84),
    ("Disulfide-linked motif", 71, 72),
    ("Leucine zipper domain",  85, 110),
]

EXPECTED_REF_LEN = 110


def read_fasta(path):
    records, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].strip(), []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def analyse(path, ref_key):
    """Return (per_position, metadata) for one aligned FASTA."""
    records = read_fasta(path)
    if not records:
        raise ValueError(f"No sequences parsed from {path}")

    seqs = [s.upper() for _, s in records]
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(
            f"{Path(path).name} is not an alignment: record lengths differ "
            f"{sorted(lengths)}"
        )
    ncol = lengths.pop()
    nseq = len(seqs)

    hits = [i for i, (h, _) in enumerate(records) if ref_key in h]
    if len(hits) != 1:
        raise ValueError(
            f"Reference '{ref_key}' matched {len(hits)} headers in "
            f"{Path(path).name}; pass an unambiguous accession with --ref."
        )
    ref_idx = hits[0]
    ref = seqs[ref_idx]
    ref_cols = [j for j in range(ncol) if ref[j] not in GAPS]

    rows = []
    for pos, j in enumerate(ref_cols, start=1):
        column = [s[j] for s in seqs]
        present = [c for c in column if c not in GAPS]
        n = len(present)
        top = Counter(present).most_common(1)[0][1]
        rows.append({
            "pos": pos,
            "aa": ref[j],
            "occupancy": n / nseq,
            "conservation": top / nseq,
            "identity": top / n,
        })

    meta = {
        "path": str(path),
        "nseq": nseq,
        "ncol": ncol,
        "ref_header": records[ref_idx][0],
        "ref_len": len(ref_cols),
    }
    return rows, meta


def domain_means(rows):
    by_pos = {r["pos"]: r for r in rows}
    out = []
    for name, start, end in DOMAINS:
        idx = [p for p in range(start, end + 1) if p in by_pos]
        if not idx:
            out.append((name, start, end, 0, None, None, None))
            continue
        k = len(idx)
        out.append((
            name, start, end, k,
            sum(by_pos[p]["occupancy"] for p in idx) / k,
            sum(by_pos[p]["conservation"] for p in idx) / k,
            sum(by_pos[p]["identity"] for p in idx) / k,
        ))
    return out


def print_report(label, note, rows, meta, per_position):
    print(f"\n{'=' * 78}")
    print(f"{label} — {note}")
    print("=" * 78)
    print(f"File             : {meta['path']}")
    print(f"Sequences        : {meta['nseq']}")
    print(f"Alignment columns: {meta['ncol']}")
    print(f"Reference        : {meta['ref_header']}")
    print(f"Reference length : {meta['ref_len']} residues", end="")
    if meta["ref_len"] != EXPECTED_REF_LEN:
        print(f"   ** expected {EXPECTED_REF_LEN} — check the alignment **")
    else:
        print()
    print()

    if per_position:
        print(f"{'Pos':>4}  {'AA':>2}  {'Occ':>7}  {'Cons':>7}  {'IdPres':>7}")
        for r in rows:
            print(f"{r['pos']:>4}  {r['aa']:>2}  {r['occupancy']:>7.3f}  "
                  f"{r['conservation']:>7.3f}  {r['identity']:>7.3f}")
        print()

    hdr = (f"{'Domain':<24}{'Residues':>10}{'Positions':>11}"
           f"{'Occupancy':>11}{'Conservation':>14}{'Id. where present':>19}")
    print(hdr)
    print("-" * len(hdr))
    for name, start, end, k, occ, con, idp in domain_means(rows):
        if occ is None:
            print(f"{name:<24}{f'{start}-{end}':>10}{'0':>11}"
                  f"{'(outside reference)':>45}")
            continue
        if k != end - start + 1:
            print(f"  ! {name}: only {k} of {end - start + 1} positions present",
                  file=sys.stderr)
        print(f"{name:<24}{f'{start}-{end}':>10}{k:>11}"
              f"{occ:>11.3f}{con:>14.3f}{idp:>19.3f}")
    print()


def write_csv(label, rows, out_dir):
    stem = label.lower().replace(" ", "_")
    pos_path = out_dir / f"{stem}_per_position.csv"
    dom_path = out_dir / f"{stem}_by_domain.csv"

    with open(pos_path, "w") as fh:
        fh.write("position,residue,occupancy,conservation_score,identity_where_present\n")
        for r in rows:
            fh.write(f"{r['pos']},{r['aa']},{r['occupancy']:.6f},"
                     f"{r['conservation']:.6f},{r['identity']:.6f}\n")

    with open(dom_path, "w") as fh:
        fh.write("domain,residues,positions,occupancy,conservation_score,"
                 "identity_where_present\n")
        for name, start, end, k, occ, con, idp in domain_means(rows):
            if occ is None:
                continue
            fh.write(f"\"{name}\",{start}-{end},{k},{occ:.6f},{con:.6f},{idp:.6f}\n")

    print(f"  wrote {pos_path.name} and {dom_path.name}")


def main():
    argv = sys.argv[1:]
    per_position = "--per-position" in argv
    want_csv = "--csv" in argv

    ref_key = DEFAULT_REF
    if "--ref" in argv:
        i = argv.index("--ref")
        if i + 1 >= len(argv):
            sys.exit("--ref needs an accession")
        ref_key = argv[i + 1]
        del argv[i:i + 2]

    positional = [a for a in argv if not a.startswith("--")]
    here = Path(__file__).resolve().parent

    if positional:
        datasets = [(f"Dataset {i}", p, "") for i, p in enumerate(positional, 1)]
    else:
        datasets = [(lbl, here / fn, note) for lbl, fn, note in DATASETS]

    failures = 0
    for label, path, note in datasets:
        path = Path(path)
        if not path.exists():
            print(f"\n[{label}] missing file: {path}", file=sys.stderr)
            failures += 1
            continue
        try:
            rows, meta = analyse(path, ref_key)
        except ValueError as exc:
            print(f"\n[{label}] {exc}", file=sys.stderr)
            failures += 1
            continue
        print_report(label, note, rows, meta, per_position)
        if want_csv:
            write_csv(label, rows, here)

    if failures:
        sys.exit(f"\n{failures} dataset(s) could not be processed.")


if __name__ == "__main__":
    main()
