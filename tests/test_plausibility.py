"""Tests for the artifact rules and the ranking that depends on them.

This is the part of the tool that decides whether a "pathogenic" call is a
finding or a genotyping error, so it is the part where being wrong does the
most damage. The thresholds are asserted as behaviour, not restated as
constants, so a change to CHIP_ERROR_RATE or COMMON_VARIANT_CEILING has to be
made deliberately rather than by accident.
"""

from __future__ import annotations

import pytest

from allele.annotate import Finding, Report, is_material
from allele.model import Sample
from allele.sources.base import (
    ZYGOSITY_ABSENT,
    ZYGOSITY_HEMI,
    ZYGOSITY_HET,
    ZYGOSITY_HOM,
    ZYGOSITY_UNKNOWN,
    Annotation,
    plausibility,
    plausibility_common_pathogenic,
    plausibility_unknown_frequency,
    stars_for,
    zygosity_for,
)


def clinvar(significance="Pathogenic", stars=2, zygosity=ZYGOSITY_HET, **kw):
    return Annotation(
        rsid="rs1", source="clinvar", category="clinical",
        significance=significance, review_stars=stars, zygosity=zygosity, **kw,
    )


# ------------------------------------------------- Hardy-Weinberg check ----

@pytest.mark.parametrize(
    "zygosity,frequency,implausible",
    [
        # Homozygous for a 1% allele: expected in 1 in 10,000, well under the
        # array's own error rate, so an error is the better explanation.
        (ZYGOSITY_HOM, 0.01, True),
        # Homozygous for a 5% allele: expected in 1 in 400. Real.
        (ZYGOSITY_HOM, 0.05, False),
        # Carrying one copy of a 1% allele is ordinary.
        (ZYGOSITY_HET, 0.01, False),
        # One copy of a 1-in-10,000 allele is not.
        (ZYGOSITY_HET, 0.0001, True),
        # On a haploid chromosome a single copy is full dosage, so the bar is f
        # itself rather than 2f(1-f).
        (ZYGOSITY_HEMI, 0.01, False),
        (ZYGOSITY_HEMI, 0.0001, True),
        # Not carrying it, or not knowing, is never an artifact claim.
        (ZYGOSITY_ABSENT, 0.0001, False),
        (ZYGOSITY_UNKNOWN, 0.0001, False),
    ],
)
def test_plausibility(zygosity, frequency, implausible):
    flagged, explanation = plausibility(zygosity, frequency)
    assert flagged is implausible
    assert (explanation is not None) is implausible


def test_plausibility_declines_to_judge_without_a_frequency():
    """No frequency is not evidence of rarity; a different rule handles that."""
    assert plausibility(ZYGOSITY_HOM, None) == (False, None)
    assert plausibility(ZYGOSITY_HOM, 0) == (False, None)


def test_plausibility_explanation_shows_its_arithmetic():
    _, explanation = plausibility(ZYGOSITY_HOM, 0.001)
    # The report's whole stance is that a demotion has to show its reasoning
    # rather than assert a verdict.
    assert "1 in 1,000,000" in explanation
    assert "%" in explanation


# ------------------------------------------------- common-pathogenic -------

def test_common_pathogenic_is_flagged():
    flagged, explanation = plausibility_common_pathogenic(0.2, "Pathogenic")
    assert flagged
    assert "20.0%" in explanation


@pytest.mark.parametrize(
    "frequency,significance",
    [
        (0.05, "Pathogenic"),      # exactly at the ceiling, not over it
        (0.2, "Benign"),           # common and benign is just a common variant
        (None, "Pathogenic"),      # nothing to judge
        (0.2, None),
    ],
)
def test_common_pathogenic_stays_quiet(frequency, significance):
    assert plausibility_common_pathogenic(frequency, significance) == (False, None)


# ------------------------------------------------- cohort-absence ----------

@pytest.mark.parametrize(
    "zygosity,significance,flagged",
    [
        (ZYGOSITY_HOM, "Pathogenic", True),
        (ZYGOSITY_HEMI, "Likely pathogenic", True),
        # A heterozygous carrier of a genuinely rare pathogenic variant is the
        # expected way to see one, so the rule must not fire here.
        (ZYGOSITY_HET, "Pathogenic", False),
        (ZYGOSITY_HOM, "Benign", False),
        (ZYGOSITY_HOM, None, False),
    ],
)
def test_unknown_frequency(zygosity, significance, flagged):
    fired, explanation = plausibility_unknown_frequency(zygosity, significance)
    assert fired is flagged
    if flagged:
        assert "1 in 100,000" in explanation


# ------------------------------------------------- zygosity ----------------

@pytest.mark.parametrize(
    "genotype,allele,expected",
    [
        ("AA", "A", ZYGOSITY_HOM),
        ("AG", "A", ZYGOSITY_HET),
        ("GG", "A", ZYGOSITY_ABSENT),
        ("A", "A", ZYGOSITY_HEMI),     # collapsed Y/MT call
        ("", "A", ZYGOSITY_UNKNOWN),
        ("AA", "", ZYGOSITY_UNKNOWN),
        ("AA", "AT", ZYGOSITY_UNKNOWN),  # multi-base alleles are not comparable
    ],
)
def test_zygosity_for(genotype, allele, expected):
    assert zygosity_for(genotype, allele) == expected


# ------------------------------------------------- review stars ------------

@pytest.mark.parametrize(
    "status,stars",
    [
        ("practice_guideline", 4),
        ("reviewed_by_expert_panel", 3),
        ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
        ("criteria_provided,_single_submitter", 1),
        ("no_assertion_criteria_provided", 0),
        # ClinVar writes these with spaces in some releases.
        ("reviewed by expert panel", 3),
        (None, None),
    ],
)
def test_stars_for(status, stars):
    assert stars_for(status) == stars


def test_unrecognized_review_status_is_unknown_not_zero():
    """A status we cannot map must not silently read as a 0-star claim."""
    assert stars_for("some_future_clinvar_status") is None


# ------------------------------------------------- ranking -----------------

def test_combined_pathogenic_classification_outranks_a_risk_factor():
    """ClinVar's combined call is a pathogenic-family result, not an unknown.

    It was previously absent from SEVERITY, so it fell to the unrecognized
    default and sorted below "risk factor" despite being the more serious
    finding.
    """
    combined = Finding("rs1", "AG", annotations=[clinvar("Pathogenic/Likely pathogenic")])
    risk = Finding("rs2", "AG", annotations=[clinvar("risk factor")])
    assert combined.score() > risk.score()


def test_severity_order_is_monotonic():
    scores = [
        Finding("rs", "AG", annotations=[clinvar(sig)]).score()
        for sig in ("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic")
    ]
    assert scores == sorted(scores, reverse=True)


def test_review_stars_scale_the_score():
    weak = Finding("rs1", "AG", annotations=[clinvar(stars=1)]).score()
    strong = Finding("rs2", "AG", annotations=[clinvar(stars=4)]).score()
    assert strong > weak


def test_an_uncarried_clinvar_entry_scores_nothing():
    absent = Finding("rs1", "GG", annotations=[clinvar(zygosity=ZYGOSITY_ABSENT)])
    assert absent.score() == 0


# ------------------------------------------------- materiality -------------

def test_zero_star_clinvar_alone_is_not_material():
    """One unreviewed submitter with no criteria is not worth a line."""
    assert not is_material(Finding("rs1", "AG", annotations=[clinvar(stars=0)]))


def test_one_star_clinvar_is_material():
    assert is_material(Finding("rs1", "AG", annotations=[clinvar(stars=1)]))


def test_snpedia_needs_magnitude():
    boring = Annotation(rsid="rs1", source="snpedia", category="curated",
                        magnitude=0.0, zygosity=ZYGOSITY_HET)
    interesting = Annotation(rsid="rs1", source="snpedia", category="curated",
                             magnitude=2.5, zygosity=ZYGOSITY_HET)
    assert not is_material(Finding("rs1", "AG", annotations=[boring]))
    assert is_material(Finding("rs1", "AG", annotations=[interesting]))


def test_uncarried_annotations_never_make_a_finding_material():
    assert not is_material(
        Finding("rs1", "GG", annotations=[clinvar(stars=4, zygosity=ZYGOSITY_ABSENT)])
    )


# ------------------------------------------------- report accounting -------

def test_headline_counts_exclude_demoted_artifacts():
    """The summary must agree with what the report displays.

    A "clinical" count that includes variants shown under the artifact heading,
    without saying so, is the number most likely to alarm someone.
    """
    real = Finding("rs1", "AG", annotations=[clinvar()])
    artifact = Finding(
        "rs2", "AA",
        annotations=[clinvar(zygosity=ZYGOSITY_HOM, flags=("implausible: see arithmetic",))],
    )
    report = Report(sample=Sample(), findings=[real, artifact])

    assert artifact.implausible and not real.implausible
    assert report.summary()["clinical"] == 1
    assert report.summary()["likely_artifacts"] == 1
    assert len(report.by_category("clinical", credible_only=False)) == 2
