"""GWAS Catalog: published trait associations, from EBI/NHGRI. Free, attribution.

This is association data, not clinical interpretation. An entry means a study
found a statistical link between a variant and a trait in some population; it
does not mean the variant causes the trait in you. Effect sizes are typically
tiny, and most associations were discovered in European-ancestry cohorts and
transfer poorly elsewhere. The report has to say so.
"""

from __future__ import annotations

import csv
import io
import re
import sqlite3
import urllib.request
import zipfile
from pathlib import Path
from typing import Iterator

from ..model import Call
from ..store import record_provenance, replace_table
from .base import Annotation, SourceInfo, zygosity_for

URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-associations_ontology-annotated-full.zip"
)
LICENSE = "GWAS Catalog (EBI/NHGRI), free with attribution"

# Genome-wide significance. Anything weaker is noise at this scale, and the
# catalog contains plenty of it.
SIGNIFICANCE_THRESHOLD = 5e-8

_COLUMNS = {
    "snps": "SNPS",
    "strongest": "STRONGEST SNP-RISK ALLELE",
    "trait": "DISEASE/TRAIT",
    "mapped": "MAPPED_TRAIT",
    "p_value": "P-VALUE",
    "effect": "OR or BETA",
    "ci": "95% CI (TEXT)",
    "genes": "MAPPED_GENE",
    "pubmed": "PUBMEDID",
    "author": "FIRST AUTHOR",
    "date": "DATE",
    "journal": "JOURNAL",
    "study": "STUDY",
    "cohort": "INITIAL SAMPLE SIZE",
    "accession": "STUDY ACCESSION",
    "risk_freq": "RISK ALLELE FREQUENCY",
}

# Ancestry labels the catalog uses in its sample descriptions. Which population
# an association was discovered in decides whether it transfers to the reader at
# all, and the great majority of the catalog is European-ancestry.
ANCESTRY_TERMS = (
    "European", "African American", "African", "East Asian", "South Asian",
    "Asian", "Hispanic", "Latino", "Native American", "Middle Eastern",
    "Oceanian", "Greater Middle Eastern", "Sub-Saharan African", "NR",
)

_SIZE_RE = re.compile(r"([0-9][0-9,]*)\s+(?:up to\s+)?[A-Za-z]")


def parse_cohort(text: str) -> tuple[int | None, tuple[str, ...]]:
    """Pull a headcount and the ancestry labels out of a sample description.

    The catalog writes these as prose, e.g. "215,551 European ancestry
    individuals, 57,332 African American individuals". Summing the numbers
    gives the discovery cohort size; the labels say who it was discovered in.
    """
    if not text:
        return None, ()
    total = 0
    for match in _SIZE_RE.finditer(text):
        try:
            total += int(match.group(1).replace(",", ""))
        except ValueError:
            pass
    found = []
    lowered = text.lower()
    for term in ANCESTRY_TERMS:
        if term.lower() in lowered and term not in found:
            # "African" would otherwise also match "African American".
            if term == "African" and "african american" in lowered and "sub-saharan" not in lowered:
                continue
            if term == "Asian" and ("east asian" in lowered or "south asian" in lowered):
                continue
            found.append(term)
    return (total or None), tuple(found)


def _risk_allele(strongest: str) -> str | None:
    """"rs4988235-T" -> "T". A "?" means the study never reported one."""
    if not strongest or "-" not in strongest:
        return None
    allele = strongest.rsplit("-", 1)[1].strip().upper()
    return allele if allele in ("A", "C", "G", "T") else None


def _iter_records(path: Path) -> Iterator[tuple]:
    with zipfile.ZipFile(path) as archive:
        members = [n for n in archive.namelist() if n.endswith((".tsv", ".txt"))]
        if not members:
            raise ValueError("GWAS Catalog archive contains no TSV")
        with archive.open(members[0]) as raw:
            stream = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            reader = csv.DictReader(stream, delimiter="\t")
            missing = [c for c in _COLUMNS.values() if c not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(f"GWAS Catalog layout changed; missing columns: {missing}")

            for row in reader:
                raw_snps = (row.get(_COLUMNS["snps"]) or "").strip()
                if not raw_snps:
                    continue

                try:
                    p_value = float(row.get(_COLUMNS["p_value"]) or "")
                except ValueError:
                    continue
                if not p_value or p_value > SIGNIFICANCE_THRESHOLD:
                    continue

                risk = _risk_allele(row.get(_COLUMNS["strongest"]) or "")
                effect = (row.get(_COLUMNS["effect"]) or "").strip()
                ci = (row.get(_COLUMNS["ci"]) or "").strip()
                if effect and ci:
                    effect = f"{effect} {ci}"

                cohort_text = (row.get(_COLUMNS["cohort"]) or "").strip()
                cohort_size, ancestries = parse_cohort(cohort_text)
                year = (row.get(_COLUMNS["date"]) or "")[:4]
                try:
                    risk_freq = float(row.get(_COLUMNS["risk_freq"]) or "")
                except ValueError:
                    risk_freq = None

                # A row can list several SNPs for one association.
                for token in raw_snps.replace(";", ",").replace(" x ", ",").split(","):
                    rsid = token.strip()
                    if not rsid.lower().startswith("rs"):
                        continue
                    yield (
                        rsid,
                        risk,
                        (row.get(_COLUMNS["trait"]) or "").strip(),
                        (row.get(_COLUMNS["mapped"]) or "").strip(),
                        p_value,
                        effect or None,
                        (row.get(_COLUMNS["genes"]) or "").strip(),
                        (row.get(_COLUMNS["pubmed"]) or "").strip(),
                        (row.get(_COLUMNS["author"]) or "").strip(),
                        year,
                        (row.get(_COLUMNS["journal"]) or "").strip(),
                        (row.get(_COLUMNS["study"]) or "").strip(),
                        (row.get(_COLUMNS["accession"]) or "").strip(),
                        cohort_text,
                        cohort_size,
                        "|".join(ancestries),
                        risk_freq,
                    )


def download(directory: Path, progress=None) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "gwas.zip"
    if progress:
        progress("downloading GWAS Catalog")
    urllib.request.urlretrieve(URL, target)
    return target


def build_index(conn: sqlite3.Connection, directory: Path, progress=None) -> int:
    path = directory / "gwas.zip"
    if not path.is_file():
        path = download(directory, progress)
    if progress:
        progress("indexing GWAS Catalog")
    columns = [
        "rsid", "risk_allele", "trait", "mapped_trait",
        "p_value", "effect_size", "genes", "pubmed_id",
        "first_author", "year", "journal", "study_title",
        "accession", "cohort", "cohort_size", "ancestries", "risk_freq",
    ]
    count = replace_table(conn, "gwas", columns, _iter_records(path))
    record_provenance(conn, "gwas", "latest", URL, LICENSE, count)
    return count


class GwasCatalog:
    """Look up published trait associations for a call."""

    name = "gwas"

    def __init__(self, conn: sqlite3.Connection, min_p: float = SIGNIFICANCE_THRESHOLD):
        self.conn = conn
        self.min_p = min_p

    def info(self) -> SourceInfo:
        row = self.conn.execute("SELECT * FROM provenance WHERE source='gwas'").fetchone()
        if not row:
            return SourceInfo(name="gwas", license=LICENSE)
        return SourceInfo(
            name="gwas",
            release=row["release"],
            downloaded=row["downloaded"],
            url=row["url"],
            license=row["license"],
            record_count=row["record_count"],
            notes=[
                "association data, not causation",
                "most studies are European-ancestry and transfer poorly",
            ],
        )

    def lookup(self, call: Call) -> list[Annotation]:
        if not call.rsid or not call.usable:
            return []
        rows = self.conn.execute(
            "SELECT * FROM gwas WHERE rsid = ? AND p_value <= ? ORDER BY p_value",
            (call.rsid, self.min_p),
        ).fetchall()

        out: list[Annotation] = []
        for row in rows:
            zygosity = zygosity_for(call.genotype, row["risk_allele"] or "", call.chrom)
            flags = ["association only, not a diagnosis"]
            if not row["risk_allele"]:
                flags.append("study reported no risk allele; carrier status unknown")

            ancestries = tuple(a for a in (row["ancestries"] or "").split("|") if a)
            if ancestries and ancestries == ("European",):
                flags.append(
                    "discovered in a European-ancestry cohort only; effect sizes "
                    "often do not transfer to other populations"
                )
            if row["cohort_size"] and row["cohort_size"] < 5000:
                flags.append(
                    f"small discovery cohort ({row['cohort_size']:,}); "
                    "small studies overstate effect sizes"
                )

            citation = None
            if row["first_author"]:
                citation = row["first_author"]
                if row["journal"]:
                    citation += f", {row['journal']}"
                if row["year"]:
                    citation += f" {row['year']}"

            out.append(
                Annotation(
                    rsid=call.rsid,
                    source="gwas",
                    category="trait",
                    title=row["trait"] or row["mapped_trait"] or "",
                    conditions=tuple(
                        t.strip() for t in (row["mapped_trait"] or "").split(",") if t.strip()
                    ),
                    genes=tuple(
                        g.strip() for g in (row["genes"] or "").split(",") if g.strip()
                    ),
                    risk_allele=row["risk_allele"],
                    p_value=row["p_value"],
                    effect_size=row["effect_size"],
                    risk_frequency=row["risk_freq"],
                    pubmed_id=row["pubmed_id"],
                    citation=citation,
                    study_title=row["study_title"] or None,
                    accession=row["accession"] or None,
                    cohort=row["cohort"] or None,
                    cohort_size=row["cohort_size"],
                    ancestries=ancestries,
                    zygosity=zygosity,
                    genotype=call.genotype,
                    flags=tuple(flags),
                )
            )
        return out
