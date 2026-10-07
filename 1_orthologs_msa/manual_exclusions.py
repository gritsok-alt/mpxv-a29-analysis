"""Manual curation exclusions and renames for Phase 5.

Edit this file — not the notebook — to curate the ortholog set. Two
dictionaries are applied in Phase 5, in this order: EXCLUSIONS first
(remove sequences), then RENAMES (relabel what remains).

EXCLUSIONS — remove accessions that pass the automated Phase 4
classification but are judged, by inspection, not to be true A29L
orthologs. Each entry MUST carry a documented reason: an accession
removed from a dataset without a recorded justification is hard to
defend in review.

    EXCLUSIONS = {
        "ACCESSION.1": "reason the sequence is excluded",
    }

RENAMES — relabel an accession to a different identifier. Use this
sparingly and only for genuine identifier corrections. Note that
malformed UniProt-style headers (sp|X|Y, tr|X|Y, db|X|Y) are cleaned
to their bare accession AUTOMATICALLY and do NOT need an entry here.

RENAMES is for the harder case: substituting one database record's
accession for another (e.g. a RefSeq NP_ accession for its GenBank
source). This CHANGES the stated provenance of the sequence — the
header will claim the new accession while the sequence bytes remain
those originally retrieved. Only do this once you have VERIFIED the
two records carry an identical sequence, and record why:

    RENAMES = {
        "OLD_ACCESSION.1": "NEW_ACCESSION.1",
    }

Both dicts may be empty ({}) if nothing is needed.
"""

EXCLUSIONS: dict[str, str] = {
    "QGT49194.1": "<Crocodylidpoxvirus below threshold GA in PF02346 QC check>",
    "QGT49408.1": "<Crocodylidpoxvirus below threshold GA in PF02346 QC check>",
    "XSG12557.1": "<Excluded due to frame-shifted C-term>",
    "QJX15523.1": "<Excluded due to being non-canonical parapoxvirus representative>",
    "VTU03264.1": "<Excluded due to wrong length, Orf virus>",
    "Q83970": "<Excluded due to wrong length, Cowpox virus>",
    "Q8QMS8": "<Excluded due to wrong length, Cowpox virus>",
    "YCL31741.1": "<Excluded due to wrong length, Cowpox virus>",
    "YFT67137.1": "<Excluded due to missing 3 first residues, Monkeypox virus>",
    "P26312": "<Excluded due to to wrong length, Vaccinia virus>",
    "YCL81897.1": "<Excluded due to 10 ambiguous residues X, Goatpox virus>",
}

RENAMES: dict[str, str] = {}
