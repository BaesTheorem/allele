"""What an annotation source is, and what one annotation looks like.

Three sources feed the report and they disagree in useful ways, so every
annotation carries its provenance and the release it came from. Nothing here
resolves a disagreement; `annotate.py` surfaces it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

from ..model import Call

# ClinVar's own review scale. The star count is the single best signal for how
# much to trust a classification: a "Pathogenic" asserted by one submitter with
# no criteria is not the same claim as one from an expert panel.
REVIEW_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_interpretation_for_the_single_variant": 0,
    "no_assertion_provided": 0,
}

# How the user's genotype relates to the allele an annotation is about.
ZYGOSITY_ABSENT = "absent"        # you do not carry the allele
ZYGOSITY_HET = "heterozygous"     # one copy
ZYGOSITY_HOM = "homozygous"       # two copies
ZYGOSITY_HEMI = "hemizygous"      # one copy on a haploid chromosome
ZYGOSITY_UNKNOWN = "unknown"      # no call, or nothing to compare against


@dataclass(frozen=True, slots=True)
class Annotation:
    """One source's statement about one variant, and how it applies to you."""

    rsid: str
    source: str                     # "clinvar" | "gwas" | "snpedia"
    category: str                   # "clinical" | "trait" | "curated"

    title: str = ""                 # human-readable headline
    conditions: tuple[str, ...] = ()
    genes: tuple[str, ...] = ()

    # Clinical (ClinVar)
    significance: str | None = None
    review_stars: int | None = None

    # Trait association (GWAS Catalog)
    risk_allele: str | None = None
    p_value: float | None = None
    effect_size: str | None = None
    pubmed_id: str | None = None

    # Curated (SNPedia, via a Promethease report)
    repute: str | None = None
    magnitude: float | None = None

    # How it applies to this person
    zygosity: str = ZYGOSITY_UNKNOWN
    genotype: str = ""
    frequency: float | None = None  # population allele frequency, 0-1

    # Caveats that must travel with the finding, never be dropped
    flags: tuple[str, ...] = ()

    @property
    def applies(self) -> bool:
        """True when the person actually carries the allele in question."""
        return self.zygosity in (ZYGOSITY_HET, ZYGOSITY_HOM, ZYGOSITY_HEMI)


@dataclass
class SourceInfo:
    """Provenance for a source, so a report can state how stale it is."""

    name: str
    release: str | None = None      # upstream release date or version
    downloaded: str | None = None   # ISO date we fetched it
    url: str | None = None
    license: str = ""
    record_count: int = 0
    notes: list[str] = field(default_factory=list)


class AnnotationSource(Protocol):
    """Every source looks the same to the annotation engine."""

    name: str

    def info(self) -> SourceInfo:
        """Where this data came from and when."""

    def lookup(self, call: Call) -> list[Annotation]:
        """Annotations for one call. Empty list if the source knows nothing."""


# Order-of-magnitude per-call error rate for consumer genotyping arrays. Used
# only to compare against how often a genotype should occur at all; the
# conclusion is robust across any plausible value here.
CHIP_ERROR_RATE = 1e-3

FLAG_IMPLAUSIBLE = "implausible"


# Smallest allele frequency the cohorts behind ClinVar's AF fields could have
# detected. ExAC alone is ~60,000 exomes, so ~120,000 alleles; a variant seen
# once would sit near 8e-6. A pathogenic variant carrying no frequency at all
# is therefore not "unknown rarity", it is "rarer than these cohorts can see".
COHORT_DETECTION_LIMIT = 1e-5


def plausibility_unknown_frequency(
    zygosity: str, significance: str | None
) -> tuple[bool, str | None]:
    """Judge a pathogenic call that carries no population frequency at all.

    ClinVar's AF fields come from 1000 Genomes, ExAC and ESP. A variant absent
    from all three is absent because it is rarer than those cohorts can
    resolve, not because nobody looked. Treating that as "frequency unknown, so
    no opinion" is how a consumer array's homozygous miscalls end up presented
    as pathogenic findings, which is the single most harmful failure mode here.
    """
    if not significance or "pathogenic" not in significance.lower():
        return False, None
    # Homozygous on an autosome, or hemizygous on a haploid chromosome (MT, Y),
    # where a single copy is already full dosage and there is no carrier state.
    if zygosity not in (ZYGOSITY_HOM, ZYGOSITY_HEMI):
        return False, None

    if zygosity == ZYGOSITY_HOM:
        expected = COHORT_DETECTION_LIMIT ** 2
        dosage = "homozygous"
        rarity = f"homozygotes would be rarer than 1 in {int(1 / expected):,}"
    else:
        expected = COHORT_DETECTION_LIMIT
        dosage = "hemizygous"
        rarity = (
            "on a haploid chromosome a single copy is full dosage, so this "
            f"would occur in fewer than 1 in {int(1 / expected):,} people"
        )

    return True, (
        f"{dosage} for a variant classified pathogenic that appears in no "
        "population frequency cohort (1000 Genomes, ExAC, ESP). That places it "
        f"below roughly 1 in {int(1 / COHORT_DETECTION_LIMIT):,} carriers, and "
        f"{rarity}. On a consumer array a genotyping error is overwhelmingly "
        "more likely. Treat as an artifact unless confirmed by clinical-grade "
        "sequencing."
    )


# A dominant pathogenic variant cannot be carried by a large fraction of the
# population; if it were, the condition would not be rare. Seeing one is the
# signature of a reference-allele or annotation error rather than a finding.
COMMON_VARIANT_CEILING = 0.05


def plausibility_common_pathogenic(
    frequency: float | None, significance: str | None
) -> tuple[bool, str | None]:
    """Catch "pathogenic" classifications on alleles most people carry."""
    if frequency is None or not significance:
        return False, None
    if "pathogenic" not in significance.lower():
        return False, None
    if frequency <= COMMON_VARIANT_CEILING:
        return False, None
    return True, (
        f"classified pathogenic yet carried by {frequency * 100:.1f}% of people. "
        "A variant that common cannot cause a rare dominant condition, so this "
        "is a reference-allele or annotation artifact rather than a finding."
    )


def plausibility(
    zygosity: str, frequency: float | None
) -> tuple[bool, str | None]:
    """Is carrying this allele more likely to be real, or a genotyping error?

    Under Hardy-Weinberg, an allele at population frequency f shows up
    homozygous in about f-squared of people and heterozygous in about 2f. When
    that expected rate falls far below the array's own error rate, the call is
    more likely an artifact than a finding. This is why a consumer chip
    reporting a rare pathogenic variant is usually wrong, and why the result
    has to be framed that way rather than as a discovery.

    Returns (implausible, explanation).
    """
    if frequency is None or frequency <= 0:
        return False, None

    if zygosity == ZYGOSITY_HOM:
        expected = frequency * frequency
    elif zygosity in (ZYGOSITY_HET, ZYGOSITY_HEMI):
        expected = 2 * frequency * (1 - frequency) if zygosity == ZYGOSITY_HET else frequency
    else:
        return False, None

    if expected >= CHIP_ERROR_RATE:
        return False, None

    one_in = int(1 / expected) if expected > 0 else 0
    ratio = CHIP_ERROR_RATE / expected if expected else float("inf")
    return True, (
        f"{zygosity} for an allele found in {frequency * 100:.3g}% of people, "
        f"expected in roughly 1 in {one_in:,}. A genotyping error is on the order "
        f"of {ratio:,.0f}x more likely than a true call. Treat as an artifact "
        f"unless confirmed by clinical-grade sequencing."
    )


def zygosity_for(genotype: str, allele: str, chrom: str | None = None) -> str:
    """Classify how many copies of `allele` a genotype carries."""
    if not genotype or not allele:
        return ZYGOSITY_UNKNOWN
    allele = allele.upper()
    if len(allele) != 1:
        return ZYGOSITY_UNKNOWN

    count = genotype.count(allele)
    if count == 0:
        return ZYGOSITY_ABSENT
    if len(genotype) == 1:
        return ZYGOSITY_HEMI
    return ZYGOSITY_HOM if count == 2 else ZYGOSITY_HET


def stars_for(review_status: str | None) -> int | None:
    """Map a ClinVar review status string to its star rating."""
    if not review_status:
        return None
    return REVIEW_STARS.get(review_status.strip().lower().replace(" ", "_"))
