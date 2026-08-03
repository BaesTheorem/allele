"""ClinVar: clinical significance for variants, from NCBI. Public domain.

We index the per-build VCF rather than variant_summary.txt.gz: it is less than
half the size, it is already keyed by position and build, and its INFO fields
carry everything the report needs.

Only single-nucleotide variants that carry an rsID survive indexing. A consumer
chip can only ever match those, so keeping the rest would inflate the database
for lookups that can never happen.
"""

from __future__ import annotations

import gzip
import sqlite3
import urllib.request
from pathlib import Path
from typing import Iterator

from ..model import Call, NUCLEOTIDES
from ..store import record_provenance, replace_table
from .base import (
    Annotation,
    FLAG_IMPLAUSIBLE,
    SourceInfo,
    ZYGOSITY_UNKNOWN,
    plausibility,
    plausibility_common_pathogenic,
    plausibility_unknown_frequency,
    stars_for,
    zygosity_for,
)

URLS = {
    37: "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
    38: "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
}
LICENSE = "Public domain (NCBI ClinVar)"

# Classifications worth putting in front of someone. Everything else in ClinVar
# is benign, uncertain, or an artifact, and reporting it as a finding would be
# noise at best and alarming at worst.
NOTABLE = (
    "pathogenic",
    "likely_pathogenic",
    "drug_response",
    "risk_factor",
)
# Substrings that disqualify a classification even if it contains "pathogenic".
NOT_NOTABLE = ("benign", "conflicting")


def _parse_info(raw: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for item in raw.split(";"):
        key, _, value = item.partition("=")
        info[key] = value
    return info


def _clean(value: str) -> str:
    """ClinVar escapes spaces as underscores and separates lists with |."""
    return value.replace("_", " ").strip()


def _iter_records(path: Path, build: int) -> Iterator[tuple]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]

            # Chip-testable variants only: one base to one base.
            if len(ref) != 1 or len(alt) != 1:
                continue
            if ref not in NUCLEOTIDES or alt not in NUCLEOTIDES:
                continue

            info = _parse_info(fields[7])
            rs = info.get("RS")
            if not rs:
                continue

            frequency = None
            for key in ("AF_TGP", "AF_EXAC", "AF_ESP"):
                if info.get(key):
                    try:
                        frequency = float(info[key])
                        break
                    except ValueError:
                        pass

            genes = ""
            if info.get("GENEINFO"):
                genes = "|".join(
                    part.split(":")[0] for part in info["GENEINFO"].split("|") if part
                )

            yield (
                f"rs{rs}",
                chrom,
                int(pos) if pos.isdigit() else None,
                build,
                ref,
                alt,
                info.get("CLNSIG", ""),
                info.get("CLNREVSTAT", ""),
                info.get("CLNDN", ""),
                genes,
                frequency,
            )


def download(build: int, directory: Path, progress=None) -> Path:
    """Fetch the ClinVar VCF for one build. Returns the local path."""
    if build not in URLS:
        raise ValueError(f"ClinVar has no VCF for build {build}")
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f"clinvar_GRCh{build}.vcf.gz"
    if progress:
        progress(f"downloading ClinVar GRCh{build}")
    urllib.request.urlretrieve(URLS[build], target)
    return target


def build_index(
    conn: sqlite3.Connection, build: int, directory: Path, progress=None
) -> int:
    """Download and index ClinVar for one build. Returns rows indexed."""
    path = download(build, directory, progress)
    if progress:
        progress("indexing ClinVar")
    columns = [
        "rsid", "chrom", "pos", "build", "ref", "alt",
        "significance", "review", "conditions", "genes", "frequency",
    ]
    count = replace_table(conn, "clinvar", columns, _iter_records(path, build))
    record_provenance(conn, "clinvar", f"GRCh{build}", URLS[build], LICENSE, count)
    return count


def is_notable(significance: str) -> bool:
    value = (significance or "").lower()
    if any(bad in value for bad in NOT_NOTABLE):
        return False
    return any(good in value for good in NOTABLE)


class ClinVar:
    """Look up clinical significance for a call."""

    name = "clinvar"

    def __init__(self, conn: sqlite3.Connection):
        self.conn = conn

    def info(self) -> SourceInfo:
        row = self.conn.execute(
            "SELECT * FROM provenance WHERE source='clinvar'"
        ).fetchone()
        if not row:
            return SourceInfo(name="clinvar", license=LICENSE)
        return SourceInfo(
            name="clinvar",
            release=row["release"],
            downloaded=row["downloaded"],
            url=row["url"],
            license=row["license"],
            record_count=row["record_count"],
        )

    def lookup(self, call: Call) -> list[Annotation]:
        if not call.rsid or not call.usable:
            return []
        rows = self.conn.execute(
            "SELECT * FROM clinvar WHERE rsid = ?", (call.rsid,)
        ).fetchall()

        out: list[Annotation] = []
        for row in rows:
            if not is_notable(row["significance"]):
                continue

            # ClinVar's ALT is the allele the classification is about. The
            # variant existing at your rsID means nothing unless you carry it.
            zygosity = zygosity_for(call.genotype, row["alt"], call.chrom)

            flags: list[str] = []
            stars = stars_for(row["review"])
            if stars is not None and stars == 0:
                flags.append("no assertion criteria: single unreviewed submission")
            if call.build and row["build"] and call.build != row["build"]:
                flags.append(f"build mismatch: call is {call.build}, ClinVar row is {row['build']}")

            implausible, explanation = plausibility_common_pathogenic(
                row["frequency"], row["significance"]
            )
            if not implausible:
                implausible, explanation = plausibility(zygosity, row["frequency"])
            if not implausible and row["frequency"] is None:
                implausible, explanation = plausibility_unknown_frequency(
                    zygosity, row["significance"]
                )
            if implausible and explanation:
                flags.append(f"{FLAG_IMPLAUSIBLE}: {explanation}")

            conditions = tuple(
                _clean(c) for c in (row["conditions"] or "").split("|")
                if c and _clean(c).lower() not in ("not specified", "not provided")
            )
            genes = tuple(g for g in (row["genes"] or "").split("|") if g)

            out.append(
                Annotation(
                    rsid=call.rsid,
                    source="clinvar",
                    category="clinical",
                    title=_clean(row["significance"]),
                    conditions=conditions,
                    genes=genes,
                    significance=_clean(row["significance"]),
                    review_stars=stars,
                    zygosity=zygosity,
                    genotype=call.genotype,
                    frequency=row["frequency"],
                    flags=tuple(flags),
                )
            )
        return out
