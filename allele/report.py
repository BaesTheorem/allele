"""Render an annotated sample to a self-contained HTML report.

No network at view time: styles are inline, there are no external fonts or
scripts, and nothing phones home. The file can be opened from a USB stick in
ten years and still work.

The design rule throughout is that a caveat travels with its finding. Flags are
never collapsed into a footnote, and anything the plausibility check demoted is
shown under its own heading explaining why, rather than being silently dropped.
"""

from __future__ import annotations

import html
from datetime import date

from .annotate import Finding, Report

CSS = """
:root{
--bg:#F8F9FF;--surface:#FFFFFF;--low:#F2F3FA;--mid:#ECEDF4;--high:#E1E2E9;
--ink:#191C20;--muted:#43474E;--line:#C3C6CF;--accent:#415F91;--accent-c:#D6E3FF;
--warn:#8B5000;--warn-c:#FFDCC1;--bad:#B3261E;--bad-c:#FFDAD6;--good:#1B6B3A;--good-c:#B6F0C6;
}
@media (prefers-color-scheme:dark){:root{
--bg:#111318;--surface:#191C20;--low:#1D2024;--mid:#282A2F;--high:#33353A;
--ink:#E1E2E9;--muted:#C3C6CF;--line:#43474E;--accent:#AAC7FF;--accent-c:#284777;
--warn:#FFB870;--warn-c:#693C00;--bad:#FFB4AB;--bad-c:#8C1D18;--good:#9BD5AB;--good-c:#005225;
}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.55 Roboto,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
-webkit-font-smoothing:antialiased}
.wrap{max-width:920px;margin:0 auto;padding:32px 20px 80px}
h1{font-size:26px;font-weight:500;margin:0 0 4px;letter-spacing:0}
h2{font-size:18px;font-weight:500;margin:36px 0 4px;letter-spacing:0}
h2 .n{color:var(--muted);font-weight:400;font-size:14px;margin-left:8px}
p.sub{color:var(--muted);margin:0 0 24px}
.banner{border:1px solid var(--line);border-left:4px solid var(--warn);
background:var(--warn-c);color:var(--ink);padding:14px 16px;margin:0 0 28px}
.banner b{display:block;margin-bottom:4px}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(120px,1fr));
gap:1px;background:var(--line);border:1px solid var(--line);margin:0 0 8px}
.cell{background:var(--low);padding:12px 14px}
.cell .v{font-size:22px;font-weight:500;font-variant-numeric:tabular-nums}
.cell .k{color:var(--muted);font-size:12px;text-transform:uppercase;letter-spacing:.07em}
.card{border:1px solid var(--line);border-left:4px solid var(--line);
background:var(--low);padding:14px 16px;margin:0 0 10px}
.card.clinical{border-left-color:var(--bad)}
.card.trait{border-left-color:var(--accent)}
.card.curated{border-left-color:var(--good)}
.card.demoted{border-left-color:var(--warn);background:var(--mid)}
.head{display:flex;flex-wrap:wrap;align-items:center;gap:8px;margin-bottom:8px}
.rsid{font-weight:500;font-size:15px}
.rsid a{color:var(--ink);text-decoration:none;border-bottom:1px solid var(--line)}
.geno{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;
background:var(--mid);border:1px solid var(--line);padding:1px 7px}
.chip{font-size:12px;padding:2px 8px;background:var(--high);color:var(--muted)}
.chip.bad{background:var(--bad-c);color:var(--ink)}
.chip.good{background:var(--good-c);color:var(--ink)}
.chip.accent{background:var(--accent-c);color:var(--ink)}
.genes{color:var(--muted);font-size:12px;margin-left:auto}
.stmt{margin:6px 0 0;padding-left:12px;border-left:2px solid var(--line)}
.stmt .src{font-size:11px;text-transform:uppercase;letter-spacing:.07em;color:var(--muted)}
.flag{margin-top:8px;padding:8px 10px;background:var(--warn-c);color:var(--ink);
font-size:12.5px;border-left:3px solid var(--warn)}
.conflict{background:var(--accent-c);border-left-color:var(--accent)}
table.src{width:100%;border-collapse:collapse;font-size:13px;margin-top:8px}
table.src td,table.src th{border-top:1px solid var(--line);padding:7px 8px;text-align:left}
table.src th{color:var(--muted);font-weight:500;font-size:11px;
text-transform:uppercase;letter-spacing:.07em;border-top:0}
footer{margin-top:48px;padding-top:16px;border-top:1px solid var(--line);
color:var(--muted);font-size:12.5px}
@media print{body{background:#fff}.card{break-inside:avoid}}
"""


def _e(value) -> str:
    return html.escape(str(value if value is not None else ""))


def _chip(text: str, kind: str = "") -> str:
    return f'<span class="chip {kind}">{_e(text)}</span>'


def _finding_card(finding: Finding, demoted: bool = False) -> str:
    categories = {a.category for a in finding.annotations}
    kind = (
        "clinical" if "clinical" in categories
        else "curated" if "curated" in categories
        else "trait"
    )
    if demoted:
        kind += " demoted"

    chips: list[str] = []
    for annotation in finding.annotations:
        if not annotation.applies:
            continue
        if annotation.significance:
            severity = "bad" if "pathogenic" in annotation.significance.lower() else "accent"
            label = annotation.significance
            if annotation.review_stars is not None:
                label += f" · {annotation.review_stars}★"
            chips.append(_chip(label, severity))
        elif annotation.magnitude:
            chips.append(_chip(f"magnitude {annotation.magnitude:g}"))
    zygosities = {a.zygosity for a in finding.annotations if a.applies}
    for zygosity in sorted(zygosities):
        chips.append(_chip(zygosity))

    # A single well-studied SNP can carry dozens of GWAS hits. Showing them all
    # buries the point, so keep the strongest few and say how many were cut.
    MAX_STATEMENTS = 5
    statements: list[str] = []
    seen: set[tuple[str, str]] = set()
    applicable = [a for a in finding.annotations if a.applies and a.title]
    applicable.sort(key=lambda a: (a.p_value is None, a.p_value or 0))
    hidden = 0
    for annotation in applicable:
        key = (annotation.source, annotation.title)
        if key in seen:
            continue
        seen.add(key)
        if len(statements) >= MAX_STATEMENTS:
            hidden += 1
            continue
        extra = ""
        if annotation.p_value:
            extra = f" &nbsp;<span class=\"src\">p={annotation.p_value:.1e}</span>"
            if annotation.effect_size:
                extra += f" &nbsp;<span class=\"src\">{_e(annotation.effect_size)}</span>"
        conditions = ""
        if annotation.conditions:
            shown = ", ".join(annotation.conditions[:4])
            conditions = f'<div class="src">{_e(shown)}</div>'
        statements.append(
            f'<div class="stmt"><span class="src">{_e(annotation.source)}</span> '
            f"{_e(annotation.title)}{extra}{conditions}</div>"
        )

    if hidden:
        statements.append(
            f'<div class="stmt"><span class="src">+{hidden} further association'
            f'{"s" if hidden != 1 else ""} not shown</span></div>'
        )

    flags = "".join(
        f'<div class="flag{" conflict" if flag in finding.conflicts else ""}">{_e(flag)}</div>'
        for flag in finding.flags
    )
    genes = (
        f'<span class="genes">{_e(", ".join(finding.genes[:4]))}</span>'
        if finding.genes else ""
    )

    return (
        f'<div class="card {kind}"><div class="head">'
        f'<span class="rsid"><a href="https://www.ncbi.nlm.nih.gov/snp/{_e(finding.rsid)}"'
        f' target="_blank" rel="noopener">{_e(finding.rsid)}</a></span>'
        f'<span class="geno">{_e(finding.genotype)}</span>'
        f'{"".join(chips)}{genes}</div>'
        f'{"".join(statements)}{flags}</div>'
    )


def render(report: Report, title: str = "Allele report") -> str:
    summary = report.summary()
    sample = report.sample

    cells = [
        ("variants read", f"{summary['calls']:,}"),
        ("annotated", f"{summary['considered']:,}"),
        ("findings", f"{len(report.credible):,}"),
        ("clinical", f"{summary['clinical']:,}"),
        ("demoted", f"{summary['likely_artifacts']:,}"),
    ]
    grid = "".join(
        f'<div class="cell"><div class="v">{_e(v)}</div><div class="k">{_e(k)}</div></div>'
        for k, v in cells
    )

    source_rows = "".join(
        f"<tr><td>{_e(s.name)}</td><td>{_e(s.release or '')}</td>"
        f"<td>{s.record_count:,}</td><td>{_e(s.downloaded or 'read from your file')}</td>"
        f"<td>{_e(s.license)}</td></tr>"
        for s in report.sources
    )

    warnings = "".join(f'<div class="flag">{_e(w)}</div>' for w in sample.warnings)

    clinical = [f for f in report.credible if any(a.category == "clinical" for a in f.annotations)]
    others = [f for f in report.credible if f not in clinical]

    def section(heading: str, items: list[Finding], demoted: bool = False, limit: int | None = None) -> str:
        if not items:
            return ""
        shown = items[:limit] if limit else items
        more = ""
        if limit and len(items) > limit:
            more = f'<p class="sub">Showing {limit:,} of {len(items):,}.</p>'
        cards = "".join(_finding_card(f, demoted) for f in shown)
        return (
            f'<h2>{_e(heading)}<span class="n">{len(items):,}</span></h2>{more}{cards}'
        )

    demoted_note = ""
    if report.artifacts:
        demoted_note = (
            '<p class="sub">These matched a database but failed a population-genetics '
            "check, so they are almost certainly genotyping errors rather than results. "
            "They are shown because hiding them entirely would be its own kind of "
            "dishonesty.</p>"
        )

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_e(title)}</title><style>{CSS}</style></head><body><div class="wrap">
<h1>{_e(title)}</h1>
<p class="sub">{_e(sample.source_format)} &middot; build {_e(sample.build or "unknown")}
&middot; generated {date.today().isoformat()}</p>

<div class="banner"><b>This is not a medical result.</b>
Consumer genotyping arrays test a fixed set of positions, miss most of the genome,
and produce false positives at rare sites. Nothing here diagnoses anything. Before
acting on any of it, confirm with a clinical-grade test and talk to a genetic
counselor or doctor.</div>

<div class="grid">{grid}</div>
{warnings}

{section("Clinical significance", clinical)}
{section("Traits and associations", others, limit=200)}

{f'<h2>Demoted as likely artifacts<span class="n">{len(report.artifacts):,}</span></h2>{demoted_note}' + "".join(_finding_card(f, True) for f in report.artifacts) if report.artifacts else ""}

<h2>Where this came from</h2>
<table class="src"><tr><th>source</th><th>release</th><th>records</th><th>obtained</th><th>license</th></tr>
{source_rows}</table>

<footer>
Generated by Allele. Your genome was read locally and never transmitted.
ClinVar is public domain. GWAS Catalog data is used with attribution to the
NHGRI-EBI Catalog of human genome-wide association studies. SNPedia content is
CC BY-NC-SA and is read from your own report rather than redistributed.
</footer>
</div></body></html>"""
