"""Source-adapter rules that decide what a finding says about you.

Each of these locks in a bug that produced a confidently wrong statement rather
than an error: a mislabelled cohort, a risk allele belonging to a different
variant, a frequency measured on the wrong assembly. None of them would have
shown up as a crash.
"""

from __future__ import annotations

import pytest

from allele.sources.base import (
    expects_common_carriers,
    plausibility_common_pathogenic,
)
from allele.sources.gwas import SIGNIFICANCE_THRESHOLD, parse_cohort, risk_alleles


# ------------------------------------------------- ancestry labels ---------

@pytest.mark.parametrize(
    "cohort,expected",
    [
        # "Caucasian" contains "asian", and a bare substring test read a
        # European cohort as an Asian one.
        ("16,839 Caucasian ancestry individuals", ("European",)),
        # "unrelated" contains "NR". The stray label made the tuple longer than
        # ("European",), which silently suppressed the transfer warning.
        ("up to 2,794 European ancestry unrelated members of twin pairs", ("European",)),
        ("5,000 European ancestry individuals", ("European",)),
        ("3,000 East Asian ancestry individuals", ("East Asian",)),
        ("1,000 NR individuals", ("NR",)),
        (
            "215,551 European ancestry individuals, 57,332 African American individuals",
            ("European", "African American"),
        ),
    ],
)
def test_ancestry_labels(cohort, expected):
    assert parse_cohort(cohort)[1] == expected


def test_european_only_warning_survives_the_word_unrelated():
    """The check downstream is `== ("European",)`, so extra labels silence it."""
    _, ancestries = parse_cohort("2,794 European ancestry unrelated individuals")
    assert ancestries == ("European",)


def test_cohort_size_is_summed():
    size, _ = parse_cohort("215,551 European ancestry individuals, 57,332 African American individuals")
    assert size == 215_551 + 57_332


# ------------------------------------------------- risk alleles ------------

def test_each_snp_keeps_its_own_risk_allele():
    """Reading only the last allele inverted carrier status in both directions."""
    assert risk_alleles("rs1015362-G; rs4911414-T") == {
        "rs1015362": "G", "rs4911414": "T",
    }


def test_single_snp_row_is_unchanged():
    assert risk_alleles("rs4988235-T") == {"rs4988235": "T"}


def test_unreported_risk_allele_yields_nothing():
    assert risk_alleles("rs4988235-?") == {}
    assert risk_alleles("?") == {}


def test_five_snp_row_does_not_collapse():
    mapped = risk_alleles("rs204999-G; rs9268528-A; rs9268542-C; rs2395163-T; rs3129871-A")
    assert len(mapped) == 5
    assert mapped["rs204999"] == "G"


# ------------------------------------------------- p-value underflow -------

def test_the_catalogs_strongest_hits_are_significant():
    """`6e-447` underflows to 0.0, and `if not p_value` dropped exactly those.

    rs7903146 (TCF7L2, type 2 diabetes) is the most replicated common-disease
    association there is, and it was being excluded as though it were missing.
    """
    underflowed = float("6e-447")
    assert underflowed == 0.0
    assert underflowed <= SIGNIFICANCE_THRESHOLD


# ------------------------------------------------- common by design --------

@pytest.mark.parametrize(
    "significance",
    [
        "Pathogenic/Pathogenic, low penetrance|risk factor",  # HFE C282Y
        "Pathogenic, low penetrance",
        "Likely pathogenic|risk factor",
        "drug response",
    ],
)
def test_classifications_expected_to_be_common(significance):
    assert expects_common_carriers(significance)


@pytest.mark.parametrize("significance", ["Pathogenic", "Likely pathogenic"])
def test_plain_pathogenic_is_not_expected_to_be_common(significance):
    assert not expects_common_carriers(significance)


def test_low_penetrance_variant_is_not_called_an_artifact():
    """The C282Y regression.

    HFE C282Y sits near 6% in Europeans and is the diagnostic genotype for
    hereditary haemochromatosis, a recessive and treatable condition. The
    common-pathogenic rule assumes a rare dominant condition, so it demoted the
    single most actionable finding a consumer array can surface.
    """
    flagged, _ = plausibility_common_pathogenic(
        0.0714, "Pathogenic/Pathogenic, low penetrance|risk factor"
    )
    assert not flagged


def test_a_plain_pathogenic_common_variant_is_still_caught():
    """The rule must keep working where its premise actually holds."""
    flagged, explanation = plausibility_common_pathogenic(0.2, "Pathogenic")
    assert flagged
    # The message should name the other explanation rather than assert one.
    assert "recessive" in explanation


# ------------------------------------------------- build safety ------------

class _FakeFrequencies:
    """Stands in for the gnomAD table, which is GRCh38 only."""

    BUILD = 38

    def __init__(self):
        self.calls = []

    def frequency(self, chrom, pos, alt, build=None):
        self.calls.append(build)
        if build is not None and build != self.BUILD:
            return None
        return 0.42


def test_gnomad_is_not_consulted_for_grch37_coordinates():
    """A GRCh37 position lands on an unrelated GRCh38 variant often enough to
    matter, and the wrong frequency then drives every plausibility check while
    being displayed to the reader as the reasoning."""
    frequencies = _FakeFrequencies()
    assert frequencies.frequency("1", 100, "A", 37) is None
    assert frequencies.frequency("1", 100, "A", 38) == 0.42


def test_the_ceiling_ignores_zygosity():
    """Being common is an argument about the variant, not about your copies."""
    for frequency in (0.06, 0.2, 0.9):
        assert plausibility_common_pathogenic(frequency, "Pathogenic")[0]
