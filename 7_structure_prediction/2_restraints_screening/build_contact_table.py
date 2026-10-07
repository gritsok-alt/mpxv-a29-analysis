#!/usr/bin/env python3
"""
Build Table S14a from a PyMOL contact listing.

Reads the output of the find_pairs loop, annotates each residue with its
heptad position, and writes TSV plus a formatted Excel workbook.

Usage:
    python build_contact_table.py contacts.txt [-o outdir]

contacts.txt is the PyMOL output pasted verbatim, e.g.

    --- A vs B: 36 contacts ---
      A/ALA  72 - B/ILE68    4.57
      ...
"""

import argparse
import re
from collections import Counter
from pathlib import Path

# Heptad register from the sequence-based assignment (Section 3.2)
A_POS = {51, 58, 65, 72, 79, 86, 93}
D_POS = {47, 54, 61, 68, 75, 82, 89, 96}


def heptad(num):
    """Position in the heptad repeat, anchored on a = 51."""
    if num in A_POS:
        return "a"
    if num in D_POS:
        return "d"
    return "abcdefg"[(num - 51) % 7]


HEADER = re.compile(r"---\s*(\w)\s+vs\s+(\w)\s*:\s*(\d+)\s+contacts")
ROW = re.compile(r"(\w)/([A-Z]{3})\s*(\d+)\s*-\s*(\w)/([A-Z]{3})\s*(\d+)\s+([\d.]+)")


def parse(path):
    rows, iface = [], None
    for line in Path(path).read_text().splitlines():
        h = HEADER.search(line)
        if h:
            iface = f"{h.group(1)}-{h.group(2)}"
            continue
        m = ROW.search(line)
        if m and iface:
            c1, n1, r1, c2, n2, r2, d = m.groups()
            r1, r2 = int(r1), int(r2)
            rows.append({
                "interface": iface,
                "chain_1": c1, "resn_1": n1.capitalize(), "resi_1": r1,
                "heptad_1": heptad(r1),
                "chain_2": c2, "resn_2": n2.capitalize(), "resi_2": r2,
                "heptad_2": heptad(r2),
                "distance": float(d),
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("contacts")
    ap.add_argument("-o", "--out", default=".")
    args = ap.parse_args()

    out = Path(args.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    rows = parse(args.contacts)
    if not rows:
        raise SystemExit("No contacts parsed - check the input format.")

    # ---- long-format table
    cols = ["interface", "chain_1", "resn_1", "resi_1", "heptad_1",
            "chain_2", "resn_2", "resi_2", "heptad_2", "distance"]
    with open(out / "TableS14a_contacts.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["interface"], r["distance"])):
            fh.write("\t".join(f"{r[c]:.2f}" if c == "distance" else str(r[c])
                               for c in cols) + "\n")

    # ---- summary by interface
    print(f"{'interface':<12}{'contacts':>10}{'a/d pairs':>12}"
          f"{'e/g involved':>14}{'b/c/f':>8}{'min':>7}{'max':>7}")
    print("-" * 70)
    summary = []
    for iface in sorted({r["interface"] for r in rows}):
        g = [r for r in rows if r["interface"] == iface]
        ad = sum(1 for r in g if r["heptad_1"] in "ad" and r["heptad_2"] in "ad")
        eg = sum(1 for r in g if "e" in (r["heptad_1"], r["heptad_2"])
                 or "g" in (r["heptad_1"], r["heptad_2"]))
        bcf = sum(1 for r in g if r["heptad_1"] in "bcf" or r["heptad_2"] in "bcf")
        d = [r["distance"] for r in g]
        print(f"{iface:<12}{len(g):>10}{ad:>12}{eg:>14}{bcf:>8}"
              f"{min(d):>7.2f}{max(d):>7.2f}")
        summary.append({"interface": iface, "contacts": len(g),
                        "a_d_pairs": ad, "e_or_g_involved": eg,
                        "b_c_f_involved": bcf,
                        "min_distance": round(min(d), 2),
                        "max_distance": round(max(d), 2)})

    # ---- which positions participate
    print("\nResidues participating, by heptad position:")
    for iface in sorted({r["interface"] for r in rows}):
        g = [r for r in rows if r["interface"] == iface]
        pos = Counter()
        for r in g:
            pos[(r["resi_1"], r["heptad_1"])] += 1
            pos[(r["resi_2"], r["heptad_2"])] += 1
        byh = {}
        for (num, h), n in pos.items():
            byh.setdefault(h, set()).add(num)
        line = "   ".join(f"{h}: {sorted(v)}"
                          for h, v in sorted(byh.items()))
        print(f"  {iface}   {line}")

    with open(out / "TableS14a_summary.tsv", "w") as fh:
        keys = list(summary[0])
        fh.write("\t".join(keys) + "\n")
        for r in summary:
            fh.write("\t".join(str(r[k]) for k in keys) + "\n")

    # ---- Excel
    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        df["pair"] = (df["chain_1"] + "/" + df["resn_1"] + df["resi_1"].astype(str)
                      + " (" + df["heptad_1"] + ")  -  "
                      + df["chain_2"] + "/" + df["resn_2"] + df["resi_2"].astype(str)
                      + " (" + df["heptad_2"] + ")")

        wide = {}
        for iface in sorted(df["interface"].unique()):
            g = df[df["interface"] == iface].sort_values("distance")
            wide[f"{iface} pair"] = list(g["pair"])
            wide[f"{iface} dist"] = list(g["distance"])
        n = max(len(v) for v in wide.values())
        for k in wide:
            wide[k] += [""] * (n - len(wide[k]))

        with pd.ExcelWriter(out / "TableS14a.xlsx", engine="openpyxl") as xw:
            pd.DataFrame(summary).to_excel(xw, sheet_name="Summary", index=False)
            pd.DataFrame(wide).to_excel(xw, sheet_name="By interface", index=False)
            df[cols].to_excel(xw, sheet_name="All contacts", index=False)
            for ws in xw.book.worksheets:
                ws.freeze_panes = "A2"
                for col in ws.columns:
                    w = max(len(str(c.value or "")) for c in col) + 2
                    ws.column_dimensions[col[0].column_letter].width = min(w, 40)
        print(f"\nExcel: {out / 'TableS14a.xlsx'}")
    except ImportError:
        print("\n(pandas/openpyxl not installed - TSV only)")

    print(f"\nWritten to {out}/")
    print("  TableS14a_contacts.tsv, TableS14a_summary.tsv, TableS14a.xlsx")


if __name__ == "__main__":
    main()
