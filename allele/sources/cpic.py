"""CPIC: pharmacogenomic prescribing guidance. Free, no registration.

CPIC publishes peer-reviewed, drug-specific guidelines for genes with enough
evidence to change prescribing. This is the most actionable thing a consumer
array can tell you, because unlike a disease risk estimate it maps to a decision
someone might already be facing.

Scope, stated plainly: this matches at the rsID level and reports that a
guideline exists for the gene. It does not call star-allele diplotypes.
Doing that properly needs phased data and full allele definitions that consumer
arrays cannot supply, and guessing a diplotype would be worse than not offering
one. PharmGKB/ClinPGx clinical annotations would add per-variant detail, but
their bulk download requires registration, so they are not used here.
"""

from __future__ import annotations

import json
import sqlite3
import urllib.request
from typing import Iterator

from ..model import Call
from ..store import record_provenance, replace_table
from .base import Annotation, SourceInfo, ZYGOSITY_HET, ZYGOSITY_HOM, ZYGOSITY_HEMI

API = "https://api.cpicpgx.org/v1"
LICENSE = "CPIC, free to use (https://cpicpgx.org)"

# CPIC's evidence levels. A and B carry prescribing actionability; C and D are
# informational and would be noise in a consumer report.
ACTIONABLE_LEVELS = ("A", "A/B", "B", "B/C")


def _get(path: str) -> list[dict]:
    request = urllib.request.Request(
        f"{API}/{path}", headers={"Accept": "application/json", "User-Agent": "allele/0.1"}
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch(progress=None) -> tuple[list[tuple], dict]:
    """Pull the CPIC tables and flatten them into per-variant rows."""
    if progress:
        progress("fetching CPIC guidelines")

    locations = _get("sequence_location?select=dbsnpid,genesymbol,name&dbsnpid=not.is.null")
    pairs = _get("pair?select=genesymbol,drugid,cpiclevel,guidelineid")
    drugs = {d["drugid"]: d["name"] for d in _get("drug?select=drugid,name")}
    guidelines = {
        g["id"]: (g.get("name"), g.get("url")) for g in _get("guideline?select=id,name,url")
    }

    by_gene: dict[str, list[tuple]] = {}
    for pair in pairs:
        level = (pair.get("cpiclevel") or "").strip()
        if level not in ACTIONABLE_LEVELS:
            continue
        drug = drugs.get(pair.get("drugid"))
        if not drug:
            continue
        name, url = guidelines.get(pair.get("guidelineid"), (None, None))
        by_gene.setdefault(pair["genesymbol"], []).append((drug, level, name, url))

    rows: list[tuple] = []
    seen: set[tuple[str, str]] = set()
    for location in locations:
        rsid = (location.get("dbsnpid") or "").strip()
        gene = (location.get("genesymbol") or "").strip()
        if not rsid.lower().startswith("rs") or gene not in by_gene:
            continue
        for drug, level, guideline, url in by_gene[gene]:
            key = (rsid, drug)
            if key in seen:
                continue
            seen.add(key)
            rows.append((rsid, gene, location.get("name") or "", drug, level, guideline or "", url or ""))

    stats = {
        "locations": len(locations),
        "genes_with_guidance": len(by_gene),
        "variants": len({r[0] for r in rows}),
    }
    return rows, stats


def build_index(conn: sqlite3.Connection, progress=None) -> int:
    rows, stats = fetch(progress)
    if progress:
        progress(f"indexing CPIC ({stats['variants']} variants, {stats['genes_with_guidance']} genes)")
    columns = ["rsid", "gene", "variant_name", "drug", "level", "guideline", "url"]
    count = replace_table(conn, "cpic", columns, iter(rows))
    record_provenance(conn, "cpic", "api", API, LICENSE, count)
    return count


class Cpic:
    """Report CPIC prescribing guidance for genes the person has variants in."""

    name = "cpic"

    def __init__(self, conn: sqlite3.Connection):
        self.conn = conn

    def info(self) -> SourceInfo:
        row = self.conn.execute("SELECT * FROM provenance WHERE source='cpic'").fetchone()
        if not row:
            return SourceInfo(name="cpic", license=LICENSE)
        return SourceInfo(
            name="cpic",
            release=row["release"],
            downloaded=row["downloaded"],
            url=row["url"],
            license=row["license"],
            record_count=row["record_count"],
            notes=[
                "gene-level guidance; star-allele diplotypes are not called",
                "levels A and B only",
            ],
        )

    def lookup(self, call: Call) -> list[Annotation]:
        if not call.rsid or not call.usable:
            return []
        rows = self.conn.execute(
            "SELECT * FROM cpic WHERE rsid = ? ORDER BY level, drug", (call.rsid,)
        ).fetchall()
        if not rows:
            return []

        # CPIC guidance is about the gene, so carriage of any genotype at a
        # guideline variant is what makes it relevant. Zygosity is reported but
        # is not a carriage test against a specific risk allele.
        if len(call.genotype) == 1:
            zygosity = ZYGOSITY_HEMI
        elif call.is_heterozygous:
            zygosity = ZYGOSITY_HET
        else:
            zygosity = ZYGOSITY_HOM

        by_gene: dict[str, list] = {}
        for row in rows:
            by_gene.setdefault(row["gene"], []).append(row)

        out: list[Annotation] = []
        for gene, entries in by_gene.items():
            drugs = sorted({e["drug"] for e in entries})
            top = entries[0]
            shown = ", ".join(drugs[:6])
            if len(drugs) > 6:
                shown += f" and {len(drugs) - 6} more"
            out.append(
                Annotation(
                    rsid=call.rsid,
                    source="cpic",
                    category="pharmacogenomic",
                    title=f"{gene}: prescribing guidance exists for {shown}",
                    conditions=tuple(drugs),
                    genes=(gene,),
                    significance=f"CPIC level {top['level']}",
                    zygosity=zygosity,
                    genotype=call.genotype,
                    citation=top["guideline"] or None,
                    accession=top["url"] or None,
                    flags=(
                        "gene-level match: this says a guideline exists for this gene, "
                        "not that your specific genotype changes a dose. Star-allele "
                        "diplotypes cannot be called from array data.",
                    ),
                )
            )
        return out
