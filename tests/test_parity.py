"""The CLI and the browser build must apply the same rules.

The README claims both implementations follow the same pipeline. That was an
assertion with nothing enforcing it, and the divergences it let through were
real: pharmacogenomic findings ranked last in the browser and first-page in the
CLI, and a measured frequency of zero flipped a finding between "credible" and
"likely artifact" depending on which one you ran.

The comparison here is rule-level, not report-level. The CLI annotates against
the full database and the browser against a shipped subset, so the two will
legitimately produce different *findings*. What must not differ is how they
judge the same input, so identical inputs go through both engines and the
verdicts are compared.

Skips cleanly when Node is unavailable, so the suite still passes on a fresh
clone without a JavaScript runtime.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from allele.annotate import Finding
from allele.model import (
    is_strand_ambiguous,
    normalize_chromosome,
    normalize_genotype,
    normalize_ploidy,
)
from allele.sources.base import (
    Annotation,
    plausibility,
    plausibility_common_pathogenic,
    plausibility_unknown_frequency,
    zygosity_for,
)

HARNESS = Path(__file__).resolve().parent / "parity" / "rules.mjs"
BUNDLE = Path(__file__).resolve().parent.parent / "docs" / "allele.js"

pytestmark = pytest.mark.skipif(
    shutil.which("node") is None or not BUNDLE.is_file(),
    reason="needs Node and docs/allele.js to compare the two implementations",
)


# Frequencies chosen to sit either side of every threshold the rules use:
# the 5% common-variant ceiling, the 1e-3 chip error rate, and zero.
FREQUENCIES = [None, 0.0, 1e-6, 1e-4, 0.001, 0.01, 0.05, 0.06, 0.2]
ZYGOSITIES = ["homozygous", "heterozygous", "hemizygous", "absent", "unknown"]
SIGNIFICANCES = [
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
    "Pathogenic, low penetrance",
    "drug response",
    "risk factor",
    None,
]

GENOTYPES = ["AG", "GA", "(C;T)", "C/T", "C|T", "--", "00", "NN", "", "II", "A", "AA"]
CHROMOSOMES = ["1", "chr1", "23", "24", "25", "26", "M", "MT", "MITO", "XY", "X", "99"]


def _cases() -> dict:
    return {
        "plausibility": [
            {"zygosity": z, "frequency": f, "significance": s}
            for z in ZYGOSITIES
            for f in FREQUENCIES
            for s in SIGNIFICANCES
        ],
        "zygosity": [
            {"genotype": g, "allele": a}
            for g in ["AA", "AG", "GG", "A", "", "AT"]
            for a in ["A", "G", "", "AT"]
        ],
        "normalizeGenotype": GENOTYPES,
        "normalizeChromosome": CHROMOSOMES,
        "normalizePloidy": [
            {"chrom": c, "genotype": g}
            for c in ["Y", "MT", "X", "1", None]
            for g in ["GG", "G", "AG", ""]
        ],
        "score": [
            [{"source": "clinvar", "significance": s, "stars": stars,
              "zygosity": "heterozygous"}]
            for s in SIGNIFICANCES for stars in [None, 0, 1, 2, 3, 4]
        ] + [
            [{"source": "cpic", "zygosity": "heterozygous"}],
            [{"source": "gwas", "zygosity": "heterozygous"}],
            [{"source": "snpedia", "magnitude": 2.5, "zygosity": "heterozygous"}],
            [{"source": "snpedia", "magnitude": 0, "zygosity": "heterozygous"}],
        ],
        "strandAmbiguous": ["AT", "CG", "AG", "AA", "TA", "GC", ""],
    }


@pytest.fixture(scope="module")
def js() -> dict:
    cases = json.dumps(_cases())
    result = subprocess.run(
        ["node", str(HARNESS)], input=cases, capture_output=True, text=True, timeout=120,
    )
    if result.returncode != 0:
        pytest.fail(f"parity harness failed:\n{result.stderr}")
    return json.loads(result.stdout)


def _python_implausible(zygosity: str, frequency: float | None, significance: str | None) -> bool:
    """Compose the three rules the way `ClinVar.lookup` composes them."""
    implausible, _ = plausibility_common_pathogenic(frequency, significance)
    if not implausible:
        implausible, _ = plausibility(zygosity, frequency)
    if not implausible and not frequency:
        implausible, _ = plausibility_unknown_frequency(zygosity, significance)
    return implausible


def test_plausibility_verdicts_match(js):
    disagreements = []
    for case, js_verdict in zip(_cases()["plausibility"], js["plausibility"]):
        mine = _python_implausible(case["zygosity"], case["frequency"], case["significance"])
        if mine != js_verdict:
            disagreements.append((case, mine, js_verdict))
    assert not disagreements, (
        f"{len(disagreements)} plausibility disagreements, first 5: {disagreements[:5]}"
    )


def test_zygosity_matches(js):
    for case, js_value in zip(_cases()["zygosity"], js["zygosity"]):
        assert zygosity_for(case["genotype"], case["allele"]) == js_value, case


def test_genotype_normalization_matches(js):
    for raw, js_value in zip(_cases()["normalizeGenotype"], js["normalizeGenotype"]):
        assert normalize_genotype(raw) == js_value, raw


def test_chromosome_normalization_matches(js):
    for raw, js_value in zip(_cases()["normalizeChromosome"], js["normalizeChromosome"]):
        assert normalize_chromosome(raw) == js_value, raw


def test_ploidy_normalization_matches(js):
    for case, js_value in zip(_cases()["normalizePloidy"], js["normalizePloidy"]):
        assert normalize_ploidy(case["chrom"], case["genotype"]) == js_value, case


def test_strand_ambiguity_matches(js):
    for genotype, js_value in zip(_cases()["strandAmbiguous"], js["strandAmbiguous"]):
        assert is_strand_ambiguous(genotype) == js_value, genotype


def test_scoring_matches(js):
    for annotations, js_value in zip(_cases()["score"], js["score"]):
        finding = Finding(
            rsid="rs1", genotype="AG",
            annotations=[
                Annotation(
                    rsid="rs1",
                    source=a["source"],
                    category="clinical",
                    significance=a.get("significance"),
                    review_stars=a.get("stars"),
                    magnitude=a.get("magnitude"),
                    zygosity=a["zygosity"],
                )
                for a in annotations
            ],
        )
        assert finding.score() == pytest.approx(js_value), annotations
