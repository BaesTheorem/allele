"""Read genotypes and SNPedia curation out of a Promethease HTML report.

Promethease is closed source, but the report it generates is self-contained
HTML with the data separated from the code: genotypes live in base64'd,
zlib-compressed JSON passed to a `decompressString(...)` call. Decoding that
recovers not just the calls but SNPedia's curation of them -- repute,
magnitude, gene, ClinVar significance -- and, uniquely among our inputs, the
strand each SNP was reported on.

That last field matters. In a real report a third of SNPs sit on the minus
strand, so a tool that joins on rsID without checking orientation will read a
third of its genotypes backwards.

Nothing here talks to Promethease or SNPedia. It reads a file already on disk.
"""

from __future__ import annotations

import base64
import json
import re
import zlib
from pathlib import Path

from .model import (
    Call,
    Sample,
    complement,
    normalize_chromosome,
    normalize_genotype,
    normalize_ploidy,
)

_PAYLOAD_RE = re.compile(r"decompressString\('([A-Za-z0-9+/=]+)'\)")

# Promethease's ClinVar significance codes, as used in its own report UI.
CLINVAR_SIGNIFICANCE = {
    1: "Untested",
    2: "Non-pathogenic",
    3: "Probable non-pathogenic",
    4: "Probable pathogenic",
    5: "Pathogenic",
    6: "Drug response",
    7: "Histocompatibility",
    255: "Other",
}


def _decode_payloads(html: str) -> list[dict]:
    """Decode every embedded JSON chunk. Bad chunks are skipped, not fatal."""
    records: list[dict] = []
    for encoded in _PAYLOAD_RE.findall(html):
        try:
            raw = zlib.decompress(base64.b64decode(encoded))
            chunk = json.loads(raw.decode("utf-8"))
        except (ValueError, zlib.error, UnicodeDecodeError):
            continue
        if isinstance(chunk, list):
            records.extend(item for item in chunk if isinstance(item, dict))
    return records


def annotation_from_record(record: dict) -> dict:
    """Pull SNPedia's curation of one SNP into a plain dict."""
    significance_code = record.get("clinvar_1")
    annotation = {
        "summary": record.get("genosummary"),
        "detail": record.get("genobody"),
        "repute": record.get("repute"),
        "magnitude": record.get("magnitude"),
        "genes": record.get("genes") or [],
        "publications": record.get("numrefs"),
        "frequency": record.get("popfreq"),
        "gmaf": record.get("gmaf"),
        "orientation": record.get("orientation"),
        "flipped": record.get("flipped"),
        "clinvar_significance": CLINVAR_SIGNIFICANCE.get(significance_code),
        "clinvar_diseases": record.get("clinvar_diseases") or [],
        "source": "snpedia",
    }
    return {k: v for k, v in annotation.items() if v not in (None, [], "")}


def genotype_reference(path: str | Path) -> dict[str, dict[str, dict]]:
    """Extract the report's SNPedia lookup table: rsID -> genotype -> grading.

    Alongside your own calls, a report embeds how SNPedia grades *every*
    genotype at the SNPs it covers, not only the one you carry. That makes a
    report you already own a usable interpretation source for a raw export
    from the same person, or a relative.

    Licensing: this is SNPedia content, CC BY-NC-SA. Reading it out of a report
    on your own disk is fine. Redistributing it is not, so this is deliberately
    a runtime extraction and nothing here is ever cached into the package.
    """
    html = Path(path).read_text(encoding="utf-8", errors="surrogateescape")
    table: dict[str, dict[str, dict]] = {}

    for encoded in _PAYLOAD_RE.findall(html):
        try:
            chunk = json.loads(zlib.decompress(base64.b64decode(encoded)).decode("utf-8"))
        except (ValueError, zlib.error, UnicodeDecodeError):
            continue
        if not isinstance(chunk, dict):
            continue
        for rsid, genotypes in chunk.items():
            if not isinstance(genotypes, dict):
                continue
            bucket = table.setdefault(rsid, {})
            for raw_genotype, grading in genotypes.items():
                if not isinstance(grading, dict):
                    continue
                key = normalize_genotype(raw_genotype)
                if not key:
                    continue
                magnitude = grading.get("mag")
                try:
                    magnitude = float(magnitude) if magnitude is not None else None
                except (TypeError, ValueError):
                    magnitude = None
                bucket[key] = {
                    "repute": grading.get("repute"),
                    "magnitude": magnitude,
                    "source": "snpedia",
                }
    return table


def parse_report(path: str | Path) -> Sample:
    """Parse a Promethease report into calls plus its SNPedia annotations."""
    path = Path(path)
    html = path.read_text(encoding="utf-8", errors="surrogateescape")
    records = _decode_payloads(html)

    if not records:
        raise ValueError(
            f"{path.name}: no Promethease genotype payload found. "
            "Is this a complete report rather than a summary page?"
        )

    calls: list[Call] = []
    orientation: dict[str, str] = {}
    snpedia_genotypes: dict[str, str] = {}
    flipped: dict[str, bool] = {}
    annotations: dict[str, dict] = {}
    builds: set[int] = set()
    genosets = 0

    for record in records:
        rsid = record.get("rsnum")
        if not rsid:
            # Genosets (gs###) are SNPedia's multi-SNP interpretations. They
            # are conclusions, not genotypes, so they are not calls.
            genosets += 1
            continue

        build = record.get("reference")
        if isinstance(build, int):
            builds.add(build)

        chrom = normalize_chromosome(record.get("chrom"))

        # `geno` is in SNPedia's orientation, which is the minus strand for
        # about a third of SNPs. `was` is the original plus-strand chip call.
        # Everything downstream (ClinVar, GWAS Catalog, other vendors' files)
        # is plus-strand, so that is what the canonical model holds; the
        # SNPedia-oriented reading is kept separately for SNPedia lookups.
        plus = normalize_genotype(record.get("was"))
        snpedia_oriented = normalize_genotype(record.get("geno"))
        if not plus:
            plus = snpedia_oriented
            if record.get("flipped") and plus:
                plus = complement(plus)

        calls.append(
            Call(
                rsid=rsid,
                chrom=chrom,
                pos=record.get("pos") if isinstance(record.get("pos"), int) else None,
                genotype=normalize_ploidy(chrom, plus),
                build=build if isinstance(build, int) else None,
            )
        )
        if snpedia_oriented:
            snpedia_genotypes[rsid] = normalize_ploidy(chrom, snpedia_oriented)

        if record.get("orientation"):
            orientation[rsid] = record["orientation"]
        flipped[rsid] = bool(record.get("flipped"))
        annotation = annotation_from_record(record)
        if annotation:
            annotations[rsid] = annotation

    warnings: list[str] = []
    if len(builds) > 1:
        warnings.append(f"report mixes genome builds {sorted(builds)}")
    if genosets:
        warnings.append(f"{genosets} genoset interpretations present (not genotypes)")

    minus = sum(1 for value in orientation.values() if value == "minus")
    if minus:
        warnings.append(
            f"{minus} of {len(orientation)} SNPs are reported on the minus strand; "
            "orientation is preserved per SNP"
        )

    return Sample(
        calls=calls,
        source_format="promethease",
        build=builds.pop() if len(builds) == 1 else None,
        path=str(path),
        warnings=warnings,
        orientation=orientation,
        flipped=flipped,
        annotations=annotations,
        alt_orientation_genotypes=snpedia_genotypes,
    )
