"""Shared NCBI E-utilities helpers.

Canonical implementation of the fetch/retry/taxonomy/protein-search
logic used across the pipeline. Phase 1 (reference database) and
Phase 3 (orphan taxa search) both import from here, so the retry
policy, pagination, and NCBI courtesy parameters are defined exactly
once - a fix applied here propagates to every phase.

NCBI courtesy parameters
------------------------
NCBI asks programmatic clients to identify themselves with an ``email``
and ``tool`` name, and to register a free ``api_key`` for anything
beyond occasional use. Supplying a key raises the request-rate ceiling
from 3 to 10 requests/second and protects long runs from IP-level
throttling. Set them once via ``configure()`` (typically from
environment variables) before calling the search functions.
"""

from __future__ import annotations

import http.client
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from xml.etree import ElementTree

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

_CREDENTIALS: dict[str, str] = {
    "email": os.environ.get("NCBI_EMAIL", ""),
    "tool": "poxvirus_pipeline",
    "api_key": os.environ.get("NCBI_API_KEY", ""),
}


def configure(email: str | None = None, tool: str | None = None,
              api_key: str | None = None) -> None:
    if email is not None:
        _CREDENTIALS["email"] = email
    if tool is not None:
        _CREDENTIALS["tool"] = tool
    if api_key is not None:
        _CREDENTIALS["api_key"] = api_key


def rate_limit_s() -> float:
    return 0.11 if _CREDENTIALS["api_key"] else 0.34


def _courtesy_params() -> dict[str, str]:
    return {k: v for k, v in _CREDENTIALS.items() if v}


def _eutils_url(endpoint: str, params: dict, eutils_base: str) -> str:
    merged = {**params, **_courtesy_params()}
    return f"{eutils_base}/{endpoint}?" + urllib.parse.urlencode(merged)


class FetchError(RuntimeError):
    """Raised when a URL could not be retrieved after all retries."""


def fetch_url(
    url: str,
    timeout: int = 300,
    max_retries: int = 5,
    retry_backoff_s: float = 3.0,
) -> bytes:
    last_exc: Exception | None = None
    for attempt in range(1, max_retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, OSError,
                http.client.IncompleteRead) as e:
            last_exc = e
            if attempt < max_retries:
                time.sleep(retry_backoff_s * attempt)
    raise FetchError(f"Failed after {max_retries} attempts: {url}") from last_exc


def get_species_taxids(
    parent_taxid: int,
    eutils_base: str = EUTILS_BASE,
    batch_size: int = 200,
    rate_limit_s: float | None = None,
) -> list[str]:
    """Return taxids of every species-level descendant of parent_taxid."""
    delay = rate_limit_s if rate_limit_s is not None else globals()["rate_limit_s"]()
    search_url = _eutils_url("esearch.fcgi", {
        "db": "taxonomy", "term": f"txid{parent_taxid}[Subtree]",
        "retmax": "100000", "retmode": "json",
    }, eutils_base)
    all_taxids = json.loads(fetch_url(search_url))["esearchresult"]["idlist"]

    species: list[str] = []
    for i in range(0, len(all_taxids), batch_size):
        batch = all_taxids[i:i + batch_size]
        summary_url = _eutils_url("esummary.fcgi", {
            "db": "taxonomy", "id": ",".join(batch), "retmode": "json",
        }, eutils_base)
        result = json.loads(fetch_url(summary_url))["result"]
        species.extend(
            tid for tid in batch
            if (rec := result.get(tid)) and rec.get("rank") == "species"
        )
        time.sleep(delay)
    return species


def esearch_protein(
    taxid: str,
    eutils_base: str = EUTILS_BASE,
    min_len: int = 60,
    max_len: int = 300,
    page_size: int = 500,
    max_ids_per_species: int | None = 3000,
    rate_limit_s: float | None = None,
    logger=None,
) -> list[str]:
    """Return NCBI protein IDs for a taxid within a length window.

    Paginates through the result set with retstart/retmax. To keep a few
    very heavily sequenced species (e.g. vaccinia, variola) from flooding
    the pool with thousands of near-identical strain deposits and stalling
    the run, pagination stops once max_ids_per_species IDs are collected.
    Set max_ids_per_species=None to retrieve every record. When the cap is
    reached the total available count is logged, so the effect is auditable.
    """
    delay = rate_limit_s if rate_limit_s is not None else globals()["rate_limit_s"]()
    term = f"txid{taxid}[Organism:exp] AND {min_len}:{max_len}[SLEN]"
    ids: list[str] = []
    retstart = 0
    while True:
        url = _eutils_url("esearch.fcgi", {
            "db": "protein", "term": term,
            "retstart": retstart, "retmax": page_size, "retmode": "json",
        }, eutils_base)
        result = json.loads(fetch_url(url))["esearchresult"]
        page = result["idlist"]
        ids.extend(page)
        retstart += len(page)
        total = int(result["count"])
        if max_ids_per_species is not None and len(ids) >= max_ids_per_species:
            ids = ids[:max_ids_per_species]
            if logger is not None:
                logger.info(
                    "taxid %s: capped at %d IDs (%d available)",
                    taxid, max_ids_per_species, total,
                )
            return ids
        if not page or retstart >= total:
            return ids
        time.sleep(delay)


def efetch_fasta(ids: list[str], eutils_base: str = EUTILS_BASE) -> str:
    """Fetch FASTA records for a batch of NCBI protein IDs."""
    url = _eutils_url("efetch.fcgi", {
        "db": "protein", "id": ",".join(ids),
        "rettype": "fasta", "retmode": "text",
    }, eutils_base)
    return fetch_url(url).decode()


def esearch_taxid_for_organism(
    organism: str,
    eutils_base: str = EUTILS_BASE,
) -> str | None:
    """Return the NCBI taxonomy id for an organism name, or None if not found."""
    url = _eutils_url("esearch.fcgi", {
        "db": "taxonomy", "term": f"{organism}[Scientific Name]",
        "retmode": "json",
    }, eutils_base)
    idlist = json.loads(fetch_url(url))["esearchresult"]["idlist"]
    return idlist[0] if idlist else None


def fetch_taxonomy_lineage(
    taxids: list[str],
    eutils_base: str = EUTILS_BASE,
    batch_size: int = 200,
    rate_limit_s: float | None = None,
) -> dict[str, dict]:
    """Map each taxid to its taxonomy record fields.

    Returns, per taxid, a dict with keys 'scientific_name', 'rank', 'genus'
    and 'genus_taxid'. Genus is read from the LineageEx ranks; when the taxon
    itself is a genus, it is reported as its own genus.
    """
    delay = rate_limit_s if rate_limit_s is not None else globals()["rate_limit_s"]()
    out: dict[str, dict] = {}
    for i in range(0, len(taxids), batch_size):
        batch = taxids[i:i + batch_size]
        url = _eutils_url("efetch.fcgi", {
            "db": "taxonomy", "id": ",".join(batch), "retmode": "xml",
        }, eutils_base)
        root = ElementTree.fromstring(fetch_url(url).decode())
        for taxon in root.findall("Taxon"):
            tid = taxon.findtext("TaxId")
            rec = {
                "scientific_name": taxon.findtext("ScientificName") or "",
                "rank": taxon.findtext("Rank") or "",
                "genus": "",
                "genus_taxid": "",
            }
            if rec["rank"] == "genus":
                rec["genus"] = rec["scientific_name"]
                rec["genus_taxid"] = tid or ""
            for node in taxon.findall("./LineageEx/Taxon"):
                if node.findtext("Rank") == "genus":
                    rec["genus"] = node.findtext("ScientificName") or ""
                    rec["genus_taxid"] = node.findtext("TaxId") or ""
                    break
            if tid:
                out[tid] = rec
        time.sleep(delay)
    return out


def fetch_organism_by_accession(accession: str, eutils_base: str = EUTILS_BASE) -> str | None:
    """Return the organism (source) name for a protein accession via NCBI efetch.

    Reads the ``/source`` or ``/organism`` qualifier from the GenPept record.
    Returns None if the record has no organism or cannot be retrieved.
    """
    url = _eutils_url("efetch.fcgi", {
        "db": "protein", "id": accession, "rettype": "gp", "retmode": "xml",
    }, eutils_base)
    try:
        root = ElementTree.fromstring(fetch_url(url).decode())
    except (FetchError, ElementTree.ParseError):
        return None
    for qual in root.iter("GBQualifier"):
        if qual.findtext("GBQualifier_name") == "organism":
            value = qual.findtext("GBQualifier_value")
            if value:
                return value.strip()
    org = root.findtext(".//GBSeq_organism")
    return org.strip() if org else None


def fetch_organism_from_uniprot(accession: str) -> str | None:
    """Return the organism name for an accession from UniProt, or None.

    Queries the UniProtKB search endpoint for the accession and reads the
    organism scientific name. Used as a fallback when NCBI has no record.
    """
    base = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": accession,
        "fields": "organism_name",
        "format": "tsv",
        "size": "1",
    }
    url = f"{base}?" + urllib.parse.urlencode(params)
    try:
        text = fetch_url(url).decode()
    except FetchError:
        return None
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        return None
    return lines[1].strip() or None
