"""SNPedia curation, read at runtime out of a Promethease report you own.

SNPedia is a wiki of human genetic variants. Its content is CC BY-NC-SA, so it
is never downloaded, cached or redistributed by this package. The one route
that respects the licence is a Promethease report already sitting on your own
disk: it embeds SNPedia's grading of every genotype at the SNPs it covers, not
just the one you carry, which makes it a usable interpretation source for a raw
export from the same person.

Orientation is the whole game here. SNPedia reports roughly a third of SNPs on
the minus strand, and the report records which. Our canonical genotypes are
plus-strand, so they are flipped into SNPedia's frame before lookup. A/T and
C/G heterozygotes read identically on both strands and are flagged rather than
guessed.
"""

from __future__ import annotations

from pathlib import Path

from ..model import Call, complement, is_strand_ambiguous
from ..promethease import genotype_reference, parse_report
from .base import (
    Annotation,
    SourceInfo,
    ZYGOSITY_HET,
    ZYGOSITY_HOM,
    ZYGOSITY_HEMI,
    ZYGOSITY_UNKNOWN,
)

LICENSE = "SNPedia content, CC BY-NC-SA. Read from your own report; never redistributed."


class SNPedia:
    """Grade a call using the SNPedia table embedded in a Promethease report."""

    name = "snpedia"

    def __init__(self, report_path: str | Path):
        self.path = Path(report_path)
        self.table = genotype_reference(self.path)
        sample = parse_report(self.path)
        self.orientation = sample.orientation
        self.flipped = sample.flipped
        self.curation = sample.annotations
        self._generated = _report_date(self.path)

    def info(self) -> SourceInfo:
        return SourceInfo(
            name="snpedia",
            release=self._generated,
            downloaded=None,
            url=str(self.path),
            license=LICENSE,
            record_count=len(self.table),
            notes=[
                "curation frozen at the moment the report was generated",
                "magnitude is a subjective interest score, not a risk estimate",
            ],
        )

    def to_snpedia_orientation(self, call: Call) -> tuple[str, bool]:
        """Flip a plus-strand genotype into SNPedia's frame for this SNP.

        Returns (genotype, ambiguous). Ambiguous means the SNP is palindromic,
        so the flip is unverifiable by allele identity.
        """
        genotype = call.genotype
        # `flipped`, not `orientation`: the two disagree on thousands of SNPs
        # and only `flipped` describes the transform actually applied.
        if self.flipped.get(call.rsid):
            if is_strand_ambiguous(genotype):
                return genotype, True
            return complement(genotype), False
        return genotype, False

    def lookup(self, call: Call) -> list[Annotation]:
        if not call.rsid or not call.usable:
            return []
        gradings = self.table.get(call.rsid)
        curation = self.curation.get(call.rsid, {})
        if not gradings and not curation:
            return []

        genotype, ambiguous = self.to_snpedia_orientation(call)
        grading = (gradings or {}).get(genotype)

        # Nothing graded for this exact genotype, and no curation either.
        if grading is None and not curation:
            return []

        repute = (grading or {}).get("repute") or curation.get("repute")
        magnitude = (grading or {}).get("magnitude")
        if magnitude is None:
            magnitude = curation.get("magnitude")

        flags: list[str] = []
        if ambiguous:
            flags.append(
                "palindromic A/T or C/G SNP on a minus-strand entry: "
                "orientation cannot be verified from the alleles"
            )
        if grading is None:
            flags.append("no grading for this exact genotype; showing SNP-level curation only")
        if magnitude is not None and magnitude >= 4:
            flags.append("high magnitude: confirm with a clinical-grade test before acting")

        if len(call.genotype) == 1:
            zygosity = ZYGOSITY_HEMI
        elif call.is_heterozygous:
            zygosity = ZYGOSITY_HET
        elif call.genotype:
            zygosity = ZYGOSITY_HOM
        else:
            zygosity = ZYGOSITY_UNKNOWN

        return [
            Annotation(
                rsid=call.rsid,
                source="snpedia",
                category="curated",
                title=curation.get("summary") or "",
                conditions=tuple(curation.get("clinvar_diseases") or []),
                genes=tuple(curation.get("genes") or []),
                significance=curation.get("clinvar_significance"),
                repute=repute,
                magnitude=magnitude,
                zygosity=zygosity,
                genotype=call.genotype,
                frequency=_as_fraction(curation.get("gmaf")),
                flags=tuple(flags),
            )
        ]


def _as_fraction(value) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if 0 <= number <= 1 else None


def _report_date(path: Path) -> str | None:
    """Pull the generation date the report stamps into its own JavaScript."""
    import re

    try:
        with open(path, "r", encoding="utf-8", errors="surrogateescape") as handle:
            head = handle.read(2_000_000)
    except OSError:
        return None
    match = re.search(
        r"generation_date\.setUTCFullYear\((\d+)\).*?setUTCMonth\((\d+)\).*?setUTCDate\((\d+)\)",
        head,
        re.S,
    )
    if not match:
        return None
    year, month, day = int(match.group(1)), int(match.group(2)) + 1, int(match.group(3))
    return f"{year:04d}-{month:02d}-{day:02d}"
