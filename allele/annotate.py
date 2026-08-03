"""Merge every source's view of a sample into one set of findings.

Sources disagree, and that is information rather than a problem to resolve. A
variant SNPedia grades as interesting can be benign in current ClinVar; a
ClinVar entry from a single unreviewed submitter can be contradicted by an
expert panel. Nothing here picks a winner. Findings keep every source's
statement, and a disagreement is surfaced as its own flag.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .model import Sample
from .sources.base import (
    FLAG_IMPLAUSIBLE,
    Annotation,
    AnnotationSource,
    SourceInfo,
)

# Order findings by how much they warrant attention. Clinical significance from
# a well-reviewed submitter outranks a subjective interest score, which outranks
# a statistical association.
SEVERITY = {
    "pathogenic": 100,
    "likely pathogenic": 90,
    "drug response": 60,
    "risk factor": 55,
}


@dataclass
class Finding:
    """Everything known about one variant, from every source that spoke."""

    rsid: str
    genotype: str
    chrom: str | None = None
    pos: int | None = None
    annotations: list[Annotation] = field(default_factory=list)
    conflicts: list[str] = field(default_factory=list)

    @property
    def sources(self) -> set[str]:
        return {a.source for a in self.annotations}

    @property
    def carried(self) -> bool:
        """Does the person actually carry the allele any source is about?"""
        return any(a.applies for a in self.annotations)

    @property
    def genes(self) -> list[str]:
        seen: list[str] = []
        for annotation in self.annotations:
            for gene in annotation.genes:
                if gene and gene not in seen:
                    seen.append(gene)
        return seen

    @property
    def flags(self) -> list[str]:
        seen: list[str] = []
        for annotation in self.annotations:
            for flag in annotation.flags:
                if flag not in seen:
                    seen.append(flag)
        return seen + self.conflicts

    @property
    def implausible(self) -> bool:
        """True when population genetics says this call is probably an artifact.

        These are kept and shown, but never ranked among real findings: a chip
        reporting a rare pathogenic variant is usually wrong, and presenting
        that as a discovery is the most harmful thing this tool could do.
        """
        return any(flag.startswith(FLAG_IMPLAUSIBLE) for flag in self.flags)

    def score(self) -> float:
        """Ranking weight. Higher sorts first."""
        best = 0.0
        for annotation in self.annotations:
            if annotation.source == "clinvar" and annotation.applies:
                weight = SEVERITY.get((annotation.significance or "").lower(), 40)
                # Review confidence scales it: a 0-star claim is not a 3-star one.
                stars = annotation.review_stars if annotation.review_stars is not None else 1
                best = max(best, weight * (0.4 + 0.15 * stars))
            elif annotation.source == "snpedia" and annotation.magnitude:
                best = max(best, float(annotation.magnitude) * 8)
            elif annotation.source == "gwas" and annotation.applies:
                best = max(best, 10)
        return best


@dataclass
class Report:
    """The annotated result for one sample."""

    sample: Sample
    findings: list[Finding] = field(default_factory=list)
    sources: list[SourceInfo] = field(default_factory=list)
    considered: int = 0

    @property
    def carried(self) -> list[Finding]:
        return [f for f in self.findings if f.carried]

    @property
    def credible(self) -> list[Finding]:
        """Findings that survive the population-genetics plausibility check."""
        return [f for f in self.findings if not f.implausible]

    @property
    def artifacts(self) -> list[Finding]:
        return [f for f in self.findings if f.implausible]

    def by_category(self, category: str) -> list[Finding]:
        return [
            f for f in self.findings
            if any(a.category == category for a in f.annotations)
        ]

    def summary(self) -> dict[str, int]:
        return {
            "calls": len(self.sample),
            "considered": self.considered,
            "findings": len(self.findings),
            "carried": len(self.carried),
            "clinical": len(self.by_category("clinical")),
            "traits": len(self.by_category("trait")),
            "curated": len(self.by_category("curated")),
            "flagged": sum(1 for f in self.findings if f.flags),
            "credible": len(self.credible),
            "likely_artifacts": len(self.artifacts),
        }


def _detect_conflicts(annotations: list[Annotation]) -> list[str]:
    """Note where sources contradict each other, without resolving it."""
    conflicts: list[str] = []

    clinical = [a for a in annotations if a.source == "clinvar" and a.applies]
    curated = [a for a in annotations if a.source == "snpedia"]

    significances = {(a.significance or "").lower() for a in clinical if a.significance}
    if len({s for s in significances if s}) > 1:
        conflicts.append(
            "ClinVar holds more than one classification for this variant: "
            + ", ".join(sorted(significances))
        )

    for curated_entry in curated:
        if curated_entry.repute == "Good" and any(
            "pathogenic" in (a.significance or "").lower() for a in clinical
        ):
            conflicts.append(
                "SNPedia grades this Good while ClinVar classifies it pathogenic; "
                "the report's curation is a snapshot and ClinVar is current"
            )
            break

    return conflicts


# SNPedia grades nearly every SNP it knows, most of them 0 for "boring". Without
# a floor, a report is tens of thousands of entries saying nothing, which buries
# the handful that matter. These are the bars for appearing at all.
MIN_MAGNITUDE = 2.0          # SNPedia's own convention: 2+ is "interesting"
MIN_REVIEW_STARS = 1         # at least one submitter provided assertion criteria


def is_material(finding: "Finding") -> bool:
    """Is this worth a line in a report a human will actually read?"""
    for annotation in finding.annotations:
        if not annotation.applies:
            continue
        if annotation.source == "clinvar":
            stars = annotation.review_stars
            if stars is None or stars >= MIN_REVIEW_STARS:
                return True
        elif annotation.source == "snpedia":
            if annotation.magnitude and annotation.magnitude >= MIN_MAGNITUDE:
                return True
        elif annotation.source == "gwas":
            # A genome-wide significant association you actually carry.
            return True
    return False


def annotate(
    sample: Sample,
    sources: list[AnnotationSource],
    include_uncarried: bool = False,
    materiality: bool = True,
) -> Report:
    """Run every source over every usable call and merge the results.

    By default only variants the person actually carries become findings.
    Reporting a pathogenic allele you do not have as a "result" is the single
    most misleading thing a tool like this can do.
    """
    findings: list[Finding] = []
    considered = 0

    for call in sample:
        if not call.usable:
            continue
        considered += 1

        annotations: list[Annotation] = []
        for source in sources:
            try:
                annotations.extend(source.lookup(call))
            except Exception as exc:  # a broken source must not lose the run
                annotations.append(
                    Annotation(
                        rsid=call.rsid,
                        source=getattr(source, "name", "unknown"),
                        category="error",
                        title=f"lookup failed: {exc}",
                    )
                )
        if not annotations:
            continue

        finding = Finding(
            rsid=call.rsid,
            genotype=call.genotype,
            chrom=call.chrom,
            pos=call.pos,
            annotations=annotations,
            conflicts=_detect_conflicts(annotations),
        )
        if not (finding.carried or include_uncarried):
            continue
        if materiality and not is_material(finding):
            continue
        findings.append(finding)

    # Artifacts sort to the bottom regardless of how severe they look.
    findings.sort(key=lambda f: (f.implausible, -f.score(), f.rsid))

    return Report(
        sample=sample,
        findings=findings,
        sources=[s.info() for s in sources],
        considered=considered,
    )
