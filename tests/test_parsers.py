"""Parser tests.

The vendor fixtures are synthetic: real 23andMe/Ancestry exports are personal
data and do not belong in a repo. They reproduce the header text and column
layout each vendor documents, plus the awkward cases (no-calls, indel markers,
numeric chromosomes, hemizygous X/Y/MT).

The Promethease test runs against a real report if one is present, and skips
otherwise, so the suite passes on a fresh clone.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from allele import parse
from allele.model import (
    Call,
    complement,
    is_strand_ambiguous,
    normalize_chromosome,
    normalize_genotype,
)

FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"


# ---------------------------------------------------------------- model ----

@pytest.mark.parametrize(
    "raw,expected",
    [
        ("AG", "AG"),
        ("GA", "AG"),          # allele order must not matter
        ("(C;T)", "CT"),       # Promethease
        ("C/T", "CT"),         # VCF-ish
        ("C|T", "CT"),
        ("  ct ", "CT"),
        ("--", ""),            # no-call spellings
        ("00", ""),
        ("NN", ""),
        ("", ""),
        ("II", "II"),          # indel markers survive normalization
        ("DD", "DD"),
        ("XY", ""),            # not nucleotides
    ],
)
def test_normalize_genotype(raw, expected):
    assert normalize_genotype(raw) == expected


@pytest.mark.parametrize(
    "raw,expected",
    [("1", "1"), ("chr7", "7"), ("23", "X"), ("24", "Y"), ("25", "X"),
     ("26", "MT"), ("M", "MT"), ("MT", "MT"), ("99", None), ("", None)],
)
def test_normalize_chromosome(raw, expected):
    assert normalize_chromosome(raw) == expected


def test_complement_round_trips():
    assert complement("AG") == "CT"
    assert complement(complement("AG")) == "AG"


def test_strand_ambiguity():
    # A/T and C/G heterozygotes read identically on either strand.
    assert is_strand_ambiguous("AT")
    assert is_strand_ambiguous("CG")
    # Everything else can be oriented by allele identity.
    assert not is_strand_ambiguous("AG")
    assert not is_strand_ambiguous("AA")
    assert not is_strand_ambiguous("CT")


def test_normalize_ploidy_collapses_haploid_chromosomes():
    from allele.model import normalize_ploidy

    # Y and MT are haploid: AncestryDNA's doubled form must collapse.
    assert normalize_ploidy("Y", "GG") == "G"
    assert normalize_ploidy("MT", "CC") == "C"
    # Autosomes and X are untouched; collapsing X would mean inferring sex.
    assert normalize_ploidy("1", "GG") == "GG"
    assert normalize_ploidy("X", "GG") == "GG"
    # A heterozygous call on a haploid chromosome is left alone rather than
    # silently halved, because it signals a miscall or heteroplasmy.
    assert normalize_ploidy("MT", "CT") == "CT"


def test_call_properties():
    assert Call("rs1", "1", 1, "").is_no_call
    assert Call("rs1", "1", 1, "II").is_indel
    assert not Call("rs1", "1", 1, "II").usable
    assert Call("rs1", "Y", 1, "G").is_hemizygous
    assert Call("rs1", "1", 1, "AG").is_heterozygous
    assert not Call("rs1", "1", 1, "AA").is_heterozygous
    assert Call("rs1", "1", 1, "AG").usable


# -------------------------------------------------------------- vendors ----

def test_23andme():
    sample = parse(FIXTURES / "23andme.txt")
    assert sample.source_format == "23andme"
    assert sample.build == 37
    calls = sample.by_rsid()
    assert calls["rs4477212"].genotype == "AA"
    assert calls["rs3094315"].genotype == "AG"
    assert calls["rs12124819"].is_no_call
    assert calls["i5000001"].is_indel
    assert calls["rs2032658"].chrom == "Y"
    assert calls["rs2853490"].chrom == "MT"
    assert calls["rs2853490"].is_hemizygous
    assert sample.stats()["no_call"] == 1
    assert sample.stats()["indel"] == 1


def test_ancestrydna_splits_alleles_and_maps_numeric_chromosomes():
    sample = parse(FIXTURES / "ancestrydna.txt")
    assert sample.source_format == "ancestrydna"
    assert sample.build == 37
    calls = sample.by_rsid()
    # Two allele columns must join, not be read as one column.
    assert calls["rs3094315"].genotype == "AG"
    assert calls["rs12124819"].is_no_call        # "0 0"
    assert calls["rs2032658"].chrom == "Y"       # 24
    assert calls["rs2853490"].chrom == "MT"      # 26
    assert calls["rs9786139"].chrom == "X"       # 23


def test_myheritage_quoted_csv():
    sample = parse(FIXTURES / "myheritage.csv")
    assert sample.source_format == "myheritage"
    calls = sample.by_rsid()
    assert calls["rs4988235"].genotype == "CT"
    assert calls["rs12124819"].is_no_call


def test_vcf_resolves_gt_indices_to_alleles():
    sample = parse(FIXTURES / "sample.vcf")
    assert sample.source_format == "vcf"
    assert sample.build == 38
    calls = sample.by_rsid()
    assert calls["rs4477212"].genotype == "AA"   # 0/0 -> REF,REF
    assert calls["rs3094315"].genotype == "AG"   # 0/1 -> REF,ALT
    assert calls["rs4988235"].genotype == "AA"   # 1/1 -> ALT,ALT
    # Multi-base ALT is an indel; it has no chip-style genotype, so it is dropped
    # rather than silently truncated to a single base.
    assert "rs99999999" not in calls
    # ./. is a no-call and is dropped.
    assert "rs88888888" not in calls


def test_vendors_agree_on_shared_snps():
    """The same SNP must normalize identically across formats."""
    shared = {}
    for name in ("23andme.txt", "ancestrydna.txt", "myheritage.csv"):
        for rsid, call in parse(FIXTURES / name).by_rsid().items():
            shared.setdefault(rsid, set()).add(call.genotype)
    for rsid, genotypes in shared.items():
        assert len(genotypes) == 1, f"{rsid} parsed inconsistently: {genotypes}"


def test_missing_build_is_reported_not_guessed():
    """A header with no build must warn, never assume 37."""
    path = FIXTURES / "_nobuild.txt"
    path.write_text("# rsid\tchromosome\tposition\tgenotype\nrs1\t1\t100\tAG\n")
    try:
        sample = parse(path)
        assert sample.build is None
        assert any("build" in w for w in sample.warnings)
    finally:
        path.unlink()


# ---------------------------------------------------------- promethease ----

REPORT = Path(os.path.expanduser("~/Downloads/promethease.original.html"))
needs_report = pytest.mark.skipif(
    not REPORT.is_file(), reason="no Promethease report available locally"
)


@needs_report
def test_promethease_report_parses():
    sample = parse(REPORT)
    assert sample.source_format == "promethease"
    assert sample.build == 37
    assert len(sample) > 10_000
    stats = sample.stats()
    assert stats["usable"] > 0
    # Promethease drops uncalled positions, so a report has no no-calls but
    # does carry the chip's raw I/D indel markers, which cannot be matched to
    # a reference allele and so must not count as usable.
    assert stats["indel"] > 0
    assert stats["usable"] == stats["total"] - stats["indel"] - stats["no_call"]


@needs_report
def test_promethease_canonical_genotypes_are_plus_strand():
    """`geno` is SNPedia-oriented; `was` is the plus-strand chip call.

    ClinVar and the GWAS Catalog are plus-strand, so the canonical model must
    hold the plus-strand reading or a third of variants mismatch.
    """
    from allele.model import complement

    sample = parse(REPORT)
    calls = sample.by_rsid()
    checked = 0
    for rsid, oriented in sample.alt_orientation_genotypes.items():
        if not sample.flipped.get(rsid):
            continue
        plus = calls[rsid].genotype
        if not plus or not oriented or is_strand_ambiguous(plus):
            continue
        # The two views must be complements of one another.
        assert complement(plus) == oriented, (rsid, plus, oriented)
        checked += 1
        if checked > 500:
            break
    assert checked > 100


@needs_report
def test_promethease_preserves_strand():
    """A third of SNPs are minus-strand; that must survive parsing."""
    sample = parse(REPORT)
    assert sample.orientation
    minus = sum(1 for v in sample.orientation.values() if v == "minus")
    assert minus > 0
    assert any("minus strand" in w for w in sample.warnings)


@needs_report
def test_promethease_genotype_reference_table():
    from allele.promethease import genotype_reference

    table = genotype_reference(REPORT)
    assert len(table) > 1_000
    # Every genotype key must already be normalized.
    for rsid, genotypes in list(table.items())[:200]:
        for genotype in genotypes:
            assert genotype == normalize_genotype(genotype)
