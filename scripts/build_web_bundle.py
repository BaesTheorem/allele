#!/usr/bin/env python3
"""Build the annotation bundle the browser version loads.

The web tool runs entirely client side, so the databases have to ship with it.
Naive JSON of everything is far too large, so this writes a columnar format
with deduplicated string tables and integer rsIDs. GitHub Pages gzips text
responses in transit, which does the rest.

Nothing SNPedia-derived is ever written here. That content is CC BY-NC-SA and
the browser reads it, if at all, out of a Promethease report the user supplies.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from allele.sources.base import stars_for  # noqa: E402
from allele.sources.clinvar import is_notable  # noqa: E402
from allele.store import connect  # noqa: E402

OUT = Path(__file__).resolve().parent.parent / "docs" / "db"

# The clinical layer is the point of the tool, so it ships whole (minus the
# unreviewed 0-star entries, which the materiality filter drops anyway).
MIN_STARS = 1
# Traits are supporting context. Keeping every association would triple the
# bundle for diminishing returns, so keep the strongest few per variant.
MAX_GWAS_PER_RSID = 3


class Strings:
    """Deduplicate repeated text into an index table."""

    def __init__(self) -> None:
        self.items: list[str] = []
        self.index: dict[str, int] = {}

    def add(self, value: str | None) -> int:
        if not value:
            return -1
        found = self.index.get(value)
        if found is None:
            found = len(self.items)
            self.index[value] = found
            self.items.append(value)
        return found


def rs_number(rsid: str) -> int | None:
    if not rsid or not rsid.lower().startswith("rs"):
        return None
    try:
        return int(rsid[2:])
    except ValueError:
        return None


def build_clinvar(conn) -> dict:
    conditions, genes = Strings(), Strings()
    significances = Strings()
    rows: list[list] = []

    query = ("SELECT rsid, alt, significance, review, conditions, genes, frequency, build, pos "
             "FROM clinvar")
    for row in conn.execute(query):
        if not is_notable(row["significance"]):
            continue
        stars = stars_for(row["review"])
        if (stars or 0) < MIN_STARS:
            continue
        number = rs_number(row["rsid"])
        if number is None or not row["alt"]:
            continue

        primary = ""
        if row["conditions"]:
            usable = [
                c.replace("_", " ").strip()
                for c in row["conditions"].split("|")
                if c and c.replace("_", " ").strip().lower()
                not in ("not specified", "not provided")
            ]
            primary = usable[0] if usable else ""
        gene = (row["genes"] or "").split("|")[0]

        rows.append([
            number,
            row["alt"],
            significances.add(row["significance"].replace("_", " ")),
            stars if stars is not None else -1,
            round(row["frequency"], 6) if row["frequency"] is not None else -1,
            conditions.add(primary),
            genes.add(gene),
            row["build"],
            row["pos"] or 0,
        ])

    rows.sort(key=lambda r: r[0])
    return {
        "columns": ["rs", "alt", "sig", "stars", "freq", "cond", "gene", "build", "pos"],
        "significance": significances.items,
        "conditions": conditions.items,
        "genes": genes.items,
        "rows": rows,
    }


def build_gwas(conn) -> dict:
    traits, genes, citations, ancestries = Strings(), Strings(), Strings(), Strings()
    per_rsid: dict[int, list] = {}

    query = (
        "SELECT rsid, risk_allele, trait, p_value, effect_size, genes, "
        "pubmed_id, first_author, year, journal, cohort_size, ancestries "
        "FROM gwas ORDER BY p_value"
    )
    for row in conn.execute(query):
        number = rs_number(row["rsid"])
        if number is None:
            continue
        bucket = per_rsid.setdefault(number, [])
        if len(bucket) >= MAX_GWAS_PER_RSID:
            continue
        trait = (row["trait"] or "").strip()
        if not trait or any(entry[1] == traits.index.get(trait) for entry in bucket):
            continue
        citation = ""
        if row["first_author"]:
            citation = row["first_author"]
            if row["journal"]:
                citation += f", {row['journal']}"
            if row["year"]:
                citation += f" {row['year']}"

        bucket.append([
            number,
            traits.add(trait),
            row["risk_allele"] or "",
            genes.add((row["genes"] or "").split(",")[0].strip()),
            row["pubmed_id"] or "",
            citations.add(citation),
            row["cohort_size"] or 0,
            ancestries.add((row["ancestries"] or "").replace("|", ", ")),
        ])

    rows = [entry for bucket in per_rsid.values() for entry in bucket]
    rows.sort(key=lambda r: r[0])
    return {
        "columns": ["rs", "trait", "risk", "gene", "pmid", "cite", "n", "ancestry"],
        "traits": traits.items,
        "genes": genes.items,
        "citations": citations.items,
        "ancestries": ancestries.items,
        "rows": rows,
    }


def build_cpic(conn) -> dict:
    """CPIC is only ~14k rows, so the browser gets all of it."""
    genes, drugs, guidelines = Strings(), Strings(), Strings()
    rows: list[list] = []
    try:
        query = "SELECT rsid, gene, drug, level, guideline, url FROM cpic"
        cursor = conn.execute(query)
    except Exception:
        return {"columns": [], "genes": [], "drugs": [], "guidelines": [], "rows": []}
    for row in cursor:
        number = rs_number(row["rsid"])
        if number is None:
            continue
        rows.append([
            number,
            genes.add(row["gene"]),
            drugs.add(row["drug"]),
            row["level"] or "",
            guidelines.add(row["url"] or ""),
        ])
    rows.sort(key=lambda r: r[0])
    return {
        "columns": ["rs", "gene", "drug", "level", "url"],
        "genes": genes.items,
        "drugs": drugs.items,
        "guidelines": guidelines.items,
        "rows": rows,
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    conn = connect()

    print("building clinvar bundle...")
    clinvar = build_clinvar(conn)
    path = OUT / "clinvar.json"
    path.write_text(json.dumps(clinvar, separators=(",", ":")), encoding="utf-8")
    print(f"  {len(clinvar['rows']):,} rows -> {path.stat().st_size / 1e6:.1f} MB")

    print("building gwas bundle...")
    gwas = build_gwas(conn)
    path = OUT / "gwas.json"
    path.write_text(json.dumps(gwas, separators=(",", ":")), encoding="utf-8")
    print(f"  {len(gwas['rows']):,} rows -> {path.stat().st_size / 1e6:.1f} MB")

    print("building cpic bundle...")
    cpic = build_cpic(conn)
    path = OUT / "cpic.json"
    path.write_text(json.dumps(cpic, separators=(",", ":")), encoding="utf-8")
    print(f"  {len(cpic['rows']):,} rows -> {path.stat().st_size / 1e6:.1f} MB")

    manifest = {
        "cpic": {
            "rows": len(cpic["rows"]),
            "license": "CPIC, free to use",
        },
        "clinvar": {
            "rows": len(clinvar["rows"]),
            "min_review_stars": MIN_STARS,
            "license": "Public domain (NCBI ClinVar)",
        },
        "gwas": {
            "rows": len(gwas["rows"]),
            "max_per_variant": MAX_GWAS_PER_RSID,
            "license": "NHGRI-EBI GWAS Catalog, free with attribution",
        },
    }
    for name, info in (("clinvar", "clinvar"), ("gwas", "gwas"), ("cpic", "cpic")):
        row = conn.execute(
            "SELECT release, downloaded FROM provenance WHERE source=?", (info,)
        ).fetchone()
        if row:
            manifest[name]["release"] = row["release"]
            manifest[name]["downloaded"] = row["downloaded"]
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"wrote manifest to {OUT / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
