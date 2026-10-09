"""Alternative output formats: JSON, annotated VCF, and report diffing.

HTML is for reading. These are for feeding something else: a pipeline, a
spreadsheet, or a comparison against a report you generated last month after
ClinVar reclassified something.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import date

from .annotate import Finding, Report
from .model import CHROMOSOMES


def to_dict(report: Report) -> dict:
    """The whole report as plain data, provenance included."""
    return {
        "generated": date.today().isoformat(),
        "sample": {
            "format": report.sample.source_format,
            "build": report.sample.build,
            "path": report.sample.path,
            "stats": report.sample.stats(),
            "warnings": report.sample.warnings,
        },
        "summary": report.summary(),
        "sources": [asdict(s) for s in report.sources],
        "findings": [_finding_dict(f, False) for f in report.credible],
        "demoted": [_finding_dict(f, True) for f in report.artifacts],
    }


def _finding_dict(finding: Finding, demoted: bool) -> dict:
    return {
        "rsid": finding.rsid,
        "genotype": finding.genotype,
        "chrom": finding.chrom,
        "pos": finding.pos,
        "genes": finding.genes,
        "score": round(finding.score(), 2),
        "demoted": demoted,
        "flags": finding.flags,
        "annotations": [
            {k: (list(v) if isinstance(v, tuple) else v) for k, v in asdict(a).items()}
            for a in finding.annotations
            if a.applies
        ],
    }


def to_json(report: Report, indent: int | None = 2) -> str:
    return json.dumps(to_dict(report), indent=indent, default=str)


def _alleles_for(finding: Finding) -> tuple[str, str] | None:
    """The REF and ALT this finding is about, if any source knows them.

    ClinVar is preferred over a prediction, since its ALT is the allele the
    classification is actually about.
    """
    ordered = sorted(
        (a for a in finding.annotations if a.ref and a.alt),
        key=lambda a: 0 if a.source == "clinvar" else 1,
    )
    for annotation in ordered:
        ref, alt = annotation.ref.upper(), annotation.alt.upper()
        if ref in _BASES and alt in _BASES:
            return ref, alt
    return None


_BASES = frozenset("ACGTN")


def _genotype_field(genotype: str, ref: str, alt: str) -> str:
    """The array's call as VCF allele indices.

    GT is a reserved, index-valued field: 0 is REF, 1 is the first ALT. Writing
    bases there (`A/G`) parses as a malformed genotype in every VCF reader.
    """
    if not genotype:
        return "./."
    indices = []
    for base in genotype:
        if base == ref:
            indices.append("0")
        elif base == alt:
            indices.append("1")
        else:
            indices.append(".")
    # Unphased genotypes are conventionally written low index first, so a
    # heterozygote reads 0/1 rather than 1/0.
    indices.sort(key=lambda i: (i == ".", i))
    return "/".join(indices)


def _contig_order(chrom: str) -> int:
    """Sort key putting chromosomes in reference order, not lexical order."""
    try:
        return CHROMOSOMES.index(chrom)
    except ValueError:
        return len(CHROMOSOMES)


def to_vcf(report: Report) -> str:
    """Findings as an annotated VCF, for feeding downstream tools.

    Only findings with a position and a known REF/ALT are emitted: a VCF row
    without trustworthy coordinates or alleles is worse than no row, because
    downstream tools cannot tell it apart from a real one.
    """
    sample = report.sample
    rows = []
    for finding in report.findings:
        if not finding.chrom or not finding.pos:
            continue
        alleles = _alleles_for(finding)
        if alleles is None:
            continue
        rows.append((finding, *alleles))

    # VCF must be coordinate-sorted within each contig to be indexable. The
    # report itself is ordered by severity, which is the right order for a human
    # and the wrong one for a file another tool will read.
    rows.sort(key=lambda item: (_contig_order(item[0].chrom), item[0].pos))

    lines = [
        "##fileformat=VCFv4.2",
        f"##fileDate={date.today().strftime('%Y%m%d')}",
        "##source=allele",
    ]
    if sample.build:
        lines.append(f"##reference=GRCh{sample.build}")
    for info in report.sources:
        lines.append(
            f"##allele_source=<name={info.name},release={info.release or 'NA'},"
            f"records={info.record_count},license=\"{info.license}\">"
        )
    # Every contig used has to be declared before the header line.
    for chrom in sorted({f.chrom for f, _, _ in rows}, key=_contig_order):
        lines.append(f"##contig=<ID={chrom}>")
    if not rows:
        lines.append(
            '##allele_note=<msg="no findings carried usable coordinates and alleles; '
            'the input may not declare a genome build">'
        )
    lines += [
        '##FILTER=<ID=DEMOTED,Description="Failed a plausibility check; likely a genotyping artifact">',
        '##INFO=<ID=AL_SRC,Number=.,Type=String,Description="Sources reporting this variant">',
        '##INFO=<ID=AL_SIG,Number=.,Type=String,Description="Clinical significance">',
        '##INFO=<ID=AL_STARS,Number=1,Type=Integer,Description="ClinVar review stars">',
        '##INFO=<ID=AL_GENE,Number=.,Type=String,Description="Genes">',
        '##INFO=<ID=AL_TRAIT,Number=.,Type=String,Description="Trait associations">',
        '##INFO=<ID=AL_DEMOTED,Number=0,Type=Flag,Description="Failed a plausibility check; likely a genotyping artifact">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype as reported by the array">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]

    def escape(value: str) -> str:
        """Make a free-text value safe as one INFO value.

        Commas separate values within a field and semicolons separate fields,
        so both have to go, along with the `=` that splits key from value.
        """
        return (
            str(value)
            .replace(" ", "_").replace(";", "|").replace(",", "|")
            .replace("=", "-").replace("\t", "_")
        )

    for finding, ref, alt in rows:
        demoted = finding.implausible
        info: list[str] = [f"AL_SRC={','.join(sorted(finding.sources))}"]

        significances = sorted(
            {a.significance for a in finding.annotations if a.applies and a.significance}
        )
        if significances:
            info.append("AL_SIG=" + ",".join(escape(s) for s in significances))
        stars = [
            a.review_stars for a in finding.annotations
            if a.applies and a.review_stars is not None
        ]
        if stars:
            info.append(f"AL_STARS={max(stars)}")
        if finding.genes:
            info.append("AL_GENE=" + ",".join(escape(g) for g in finding.genes[:5]))
        traits = sorted(
            {a.title for a in finding.annotations if a.applies and a.category == "trait" and a.title}
        )
        if traits:
            info.append("AL_TRAIT=" + ",".join(escape(t) for t in traits[:5]))
        if demoted:
            info.append("AL_DEMOTED")

        lines.append(
            "\t".join([
                finding.chrom, str(finding.pos), finding.rsid or ".",
                ref, alt, ".", "DEMOTED" if demoted else "PASS",
                ";".join(info), "GT", _genotype_field(finding.genotype, ref, alt),
            ])
        )

    return "\n".join(lines) + "\n"


def diff(old: dict, new: dict) -> dict:
    """Compare two JSON reports for the same person.

    The point is drift in the databases, not in you: ClinVar reclassifies
    variants continuously, so a finding can appear, vanish, or change severity
    without your genome changing at all.
    """
    def index(report: dict) -> dict[str, dict]:
        out = {}
        for finding in report.get("findings", []) + report.get("demoted", []):
            out[finding["rsid"]] = finding
        return out

    before, after = index(old), index(new)
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))

    changed = []
    for rsid in sorted(set(before) & set(after)):
        was, now = before[rsid], after[rsid]
        deltas = {}
        if was.get("demoted") != now.get("demoted"):
            deltas["demoted"] = [was.get("demoted"), now.get("demoted")]
        was_sig = sorted({a.get("significance") for a in was["annotations"] if a.get("significance")})
        now_sig = sorted({a.get("significance") for a in now["annotations"] if a.get("significance")})
        if was_sig != now_sig:
            deltas["significance"] = [was_sig, now_sig]
        if deltas:
            changed.append({"rsid": rsid, "genes": now.get("genes"), "changes": deltas})

    return {
        "before": old.get("generated"),
        "after": new.get("generated"),
        "added": [after[r] for r in added],
        "removed": [before[r] for r in removed],
        "changed": changed,
        "counts": {"added": len(added), "removed": len(removed), "changed": len(changed)},
    }
