"""Normalized genotype model shared by every input format.

The whole point of this module is that a 23andMe export, an AncestryDNA
export and a Promethease report all collapse to the same `Call` records, so
annotation and reporting downstream never need to know where the data came
from.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator

# Chromosome names, normalized. Vendors disagree: AncestryDNA uses 23/24/25/26,
# 23andMe uses X/Y/MT, VCFs may prefix with "chr".
CHROMOSOMES = tuple(str(i) for i in range(1, 23)) + ("X", "Y", "MT")

_CHROM_ALIASES = {
    "23": "X",
    "24": "Y",
    "25": "X",  # AncestryDNA: pseudo-autosomal region of X
    "26": "MT",
    "M": "MT",
    "MITO": "MT",
    "XY": "X",
}

# Genotype tokens that mean "the chip did not produce a call here".
NO_CALL_TOKENS = frozenset({"--", "-", "00", "0", "NN", "N", "", "..", "."})

# Bases we can reason about. I/D are 23andMe's insertion/deletion markers and
# carry no sequence, so they can never be matched against a reference allele.
NUCLEOTIDES = frozenset("ACGT")
INDEL_MARKERS = frozenset("ID")

COMPLEMENT = str.maketrans("ACGT", "TGCA")

# A SNP whose two alleles are complements of each other (A/T or C/G) cannot be
# strand-resolved by allele identity alone: the plus-strand and minus-strand
# readings are both valid-looking. These need external orientation metadata.
AMBIGUOUS_PAIRS = frozenset({frozenset("AT"), frozenset("CG")})


def normalize_chromosome(raw: str) -> str | None:
    """Return a canonical chromosome name, or None if unrecognized."""
    if raw is None:
        return None
    value = str(raw).strip().upper()
    if value.startswith("CHR"):
        value = value[3:]
    value = _CHROM_ALIASES.get(value, value)
    return value if value in CHROMOSOMES else None


def normalize_genotype(raw: str) -> str:
    """Uppercase, strip separators, sort alleles so AG and GA compare equal.

    Returns "" for a no-call. Promethease writes "(C;T)", VCFs write "C/T" or
    "C|T", chip exports write "CT"; all three land on "CT".
    """
    if raw is None:
        return ""
    value = str(raw).strip().upper()
    for junk in "();|/, \t":
        value = value.replace(junk, "")
    if value in NO_CALL_TOKENS:
        return ""
    # Any residual non-allele character means we do not understand this call.
    if not all(ch in NUCLEOTIDES or ch in INDEL_MARKERS for ch in value):
        return ""
    return "".join(sorted(value))


# Chromosomes that are haploid wherever they are present at all.
HAPLOID_CHROMOSOMES = frozenset({"Y", "MT"})


def normalize_ploidy(chrom: str | None, genotype: str) -> str:
    """Collapse a doubled call on a haploid chromosome to one allele.

    Vendors disagree here: 23andMe writes a male Y or an MT call as "G", while
    AncestryDNA's two-allele-column layout writes the same call as "G G". Y and
    MT are haploid, so the single allele is the honest reading and the doubled
    form is a formatting artifact. Without this, the same person's Ancestry and
    23andMe files disagree on every Y and MT site.

    X is deliberately left alone: collapsing it would mean inferring sex, and
    being wrong about that is worse than the inconsistency.
    """
    if (
        chrom in HAPLOID_CHROMOSOMES
        and len(genotype) == 2
        and genotype[0] == genotype[1]
    ):
        return genotype[0]
    return genotype


def complement(genotype: str) -> str:
    """Flip a genotype to the other strand. Indel markers pass through."""
    return "".join(sorted(genotype.translate(COMPLEMENT)))


def is_strand_ambiguous(genotype: str) -> bool:
    """True for heterozygous A/T or C/G calls, which no allele check can orient."""
    alleles = set(genotype)
    return len(alleles) == 2 and frozenset(alleles) in AMBIGUOUS_PAIRS


@dataclass(frozen=True, slots=True)
class Call:
    """One genotyped position, normalized."""

    rsid: str
    chrom: str | None
    pos: int | None
    genotype: str  # "" means no-call
    build: int | None = None  # 36, 37 or 38

    @property
    def is_no_call(self) -> bool:
        return not self.genotype

    @property
    def is_indel(self) -> bool:
        return any(ch in INDEL_MARKERS for ch in self.genotype)

    @property
    def is_hemizygous(self) -> bool:
        """Single-allele call, normal on male X/Y and on mitochondrial DNA."""
        return len(self.genotype) == 1

    @property
    def is_heterozygous(self) -> bool:
        return len(self.genotype) == 2 and self.genotype[0] != self.genotype[1]

    @property
    def usable(self) -> bool:
        """Can this call be matched against a reference allele database?"""
        return bool(self.genotype) and not self.is_indel


@dataclass
class Sample:
    """A parsed input file: the calls plus everything we learned about them."""

    calls: list[Call] = field(default_factory=list)
    source_format: str = "unknown"
    build: int | None = None
    path: str | None = None
    warnings: list[str] = field(default_factory=list)
    # Populated only by sources that know the strand per SNP (Promethease does;
    # raw chip exports do not, because they are all plus-strand by convention).
    orientation: dict[str, str] = field(default_factory=dict)
    # Whether the curating source flipped this SNP relative to the chip
    # reading. This, not `orientation`, is the field that drives the
    # transform: the two disagree on thousands of SNPs.
    flipped: dict[str, bool] = field(default_factory=dict)
    # Curation that arrived with the input, keyed by rsID. Only Promethease
    # reports carry this; raw exports are annotated from local databases later.
    annotations: dict[str, dict] = field(default_factory=dict)
    # Genotypes as the curating source reported them, when that is not the plus
    # strand. `calls` is always plus-strand; this is the SNPedia-oriented view.
    alt_orientation_genotypes: dict[str, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.calls)

    def __iter__(self) -> Iterator[Call]:
        return iter(self.calls)

    def by_rsid(self) -> dict[str, Call]:
        """Last call wins. Duplicate rsIDs occur in some vendor exports."""
        return {c.rsid: c for c in self.calls if c.rsid}

    def stats(self) -> dict[str, int]:
        total = len(self.calls)
        no_call = sum(1 for c in self.calls if c.is_no_call)
        indel = sum(1 for c in self.calls if c.is_indel)
        return {
            "total": total,
            "usable": sum(1 for c in self.calls if c.usable),
            "no_call": no_call,
            "indel": indel,
            "heterozygous": sum(1 for c in self.calls if c.is_heterozygous),
            "hemizygous": sum(1 for c in self.calls if c.is_hemizygous),
        }
