"""VCF export has to be readable by the tools it exists to feed.

The previous implementation emitted `.` for REF and ALT and wrote the array's
bases into GT, which is a reserved index-valued field. Every record was rejected
by htslib, so the format was output nothing could consume. These tests assert
the spec rules that were being broken, and the pysam round-trip at the end is
the positive control: it is an independent parser, so it can fail.
"""

from __future__ import annotations

import pytest

from allele.annotate import Finding, Report
from allele.export import to_vcf
from allele.model import Sample
from allele.sources.base import Annotation


def finding(rsid, chrom, pos, genotype, ref="G", alt="A", **kw):
    return Finding(
        rsid=rsid, genotype=genotype, chrom=chrom, pos=pos,
        annotations=[Annotation(
            rsid=rsid, source="clinvar", category="clinical",
            ref=ref, alt=alt, significance="Pathogenic", review_stars=2,
            zygosity="heterozygous", **kw,
        )],
    )


def render(findings):
    return to_vcf(Report(sample=Sample(build=38), findings=findings))


def body(text):
    return [l for l in text.splitlines() if l and not l.startswith("#")]


def test_ref_and_alt_are_real_bases():
    """`.` is not a legal REF; the spec requires [ACGTN]+."""
    fields = body(render([finding("rs1", "1", 100, "GA")]))[0].split("\t")
    assert fields[3] == "G"
    assert fields[4] == "A"


def test_genotype_is_written_as_allele_indices():
    rows = body(render([
        finding("rs1", "1", 100, "GA"),   # het  -> 0/1
        finding("rs2", "1", 200, "AA"),   # hom alt -> 1/1
        finding("rs3", "1", 300, "GG"),   # hom ref -> 0/0
    ]))
    assert [r.split("\t")[9] for r in rows] == ["0/1", "1/1", "0/0"]


def test_records_are_coordinate_sorted_in_contig_order():
    """Severity order is right for a human and wrong for a file."""
    text = render([
        finding("rs1", "X", 10, "GA"),
        finding("rs2", "2", 900, "GA"),
        finding("rs3", "2", 100, "GA"),
        finding("rs4", "11", 50, "GA"),
    ])
    rows = [r.split("\t")[:2] for r in body(text)]
    # Reference order, so 2 < 11 < X rather than the lexical 11 < 2 < X.
    assert rows == [["2", "100"], ["2", "900"], ["11", "50"], ["X", "10"]]


def test_every_meta_line_precedes_the_header():
    lines = render([]).splitlines()
    header = lines.index("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")
    assert all(l.startswith("##") for l in lines[:header])


def test_contigs_are_declared():
    text = render([finding("rs1", "7", 100, "GA")])
    assert "##contig=<ID=7>" in text


def test_demoted_filter_is_declared_before_use():
    demoted = finding("rs1", "1", 100, "AA", flags=("implausible: see arithmetic",))
    text = render([demoted])
    assert "##FILTER=<ID=DEMOTED," in text
    assert body(text)[0].split("\t")[6] == "DEMOTED"


def test_findings_without_alleles_are_skipped():
    """A row downstream tools cannot interpret is worse than no row."""
    bare = Finding(
        rsid="rs1", genotype="GA", chrom="1", pos=100,
        annotations=[Annotation(rsid="rs1", source="gwas", category="trait",
                                zygosity="heterozygous")],
    )
    assert body(render([bare])) == []


def test_free_text_never_introduces_extra_info_values():
    """A significance containing a comma silently became two INFO values.

    `,` separates values inside one INFO field and `;` separates fields, so a
    single significance string carrying either has to come out as one value.
    """
    messy = finding("rs1", "1", 100, "GA")
    messy.annotations[0] = Annotation(
        rsid="rs1", source="clinvar", category="clinical", ref="G", alt="A",
        significance="Pathogenic, low penetrance; unusual", review_stars=2,
        zygosity="heterozygous",
    )
    info = dict(
        entry.partition("=")[::2] for entry in body(render([messy]))[0].split("\t")[7].split(";")
    )
    assert info["AL_SIG"].count(",") == 0, "one significance became several values"


def test_output_parses_with_an_independent_vcf_reader():
    pysam = pytest.importorskip("pysam", reason="pysam is the gnomAD extra")
    text = render([
        finding("rs1", "1", 100, "GA"),
        finding("rs2", "X", 200, "AA"),
    ])
    path = pytest.importorskip("tempfile")
    import tempfile, os
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as handle:
        handle.write(text)
        name = handle.name
    try:
        records = list(pysam.VariantFile(name))
        assert [r.ref for r in records] == ["G", "G"]
        assert [r.alts[0] for r in records] == ["A", "A"]
        assert [r.samples[0]["GT"] for r in records] == [(0, 1), (1, 1)]
    finally:
        os.unlink(name)
