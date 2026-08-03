"""Command line interface.

    allele db update              download and index ClinVar + GWAS Catalog
    allele report FILE            annotate a DNA export and write HTML
    allele info FILE              what is in this file, without annotating it
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .annotate import annotate
from .builds import verify as verify_build
from .parsers import parse
from .export import diff, to_json, to_vcf
from .report import render
from .store import DEFAULT_DIR, connect, finalize, is_populated, provenance


def _progress(message: str) -> None:
    print(f"  {message}", file=sys.stderr, flush=True)


def cmd_db(args: argparse.Namespace) -> int:
    from .sources import alphamissense, clinvar, cpic, gnomad, gwas

    directory = Path(args.data_dir) if args.data_dir else DEFAULT_DIR
    downloads = directory / "downloads"
    conn = connect(directory)

    if args.action == "path":
        print(directory)
        return 0

    if args.action == "clean":
        import shutil

        if downloads.exists():
            size = sum(f.stat().st_size for f in downloads.rglob("*") if f.is_file())
            shutil.rmtree(downloads)
            print(f"removed {size / 1e6:.0f} MB of downloaded source files from {downloads}")
        else:
            print("nothing to clean")
        print("indexes kept; run 'allele db update' to refetch sources")
        return 0

    if args.action == "status":
        rows = provenance(conn)
        if not rows:
            print("no databases indexed yet; run: allele db update")
            return 1
        for name, row in rows.items():
            print(
                f"{name:9s} release={row['release'] or '?':10s} "
                f"records={row['record_count']:>10,}  downloaded={row['downloaded']}"
            )
        return 0

    print("Updating annotation databases. This downloads a few hundred MB once.", file=sys.stderr)
    if not args.only or args.only == "clinvar":
        count = clinvar.build_index(conn, args.build, downloads, _progress)
        print(f"  ClinVar: {count:,} variants", file=sys.stderr)
    if not args.only or args.only == "gwas":
        count = gwas.build_index(conn, downloads, _progress)
        print(f"  GWAS Catalog: {count:,} associations", file=sys.stderr)
    if not args.only or args.only == "cpic":
        count = cpic.build_index(conn, _progress)
        print(f"  CPIC: {count:,} variant-drug pairs", file=sys.stderr)
    # Enrichment sources are GRCh38-only and slower, so they are opt-in.
    if args.only == "gnomad":
        count = gnomad.build_index(conn, _progress)
        print(f"  gnomAD: {count:,} frequencies", file=sys.stderr)
    if args.only == "alphamissense":
        count = alphamissense.build_index(conn, downloads, _progress)
        print(f"  AlphaMissense: {count:,} predictions", file=sys.stderr)
    finalize(conn)
    print("Done.", file=sys.stderr)
    return 0


def _build_sources(conn, snpedia_report: str | None) -> list:
    from .sources import AlphaMissense, ClinVar, Cpic, Gnomad, GwasCatalog, SNPedia

    sources: list = []
    # gnomAD is an enricher, not a source of findings: it supplies the
    # frequency the plausibility checks need rather than producing its own.
    frequencies = Gnomad(conn) if is_populated(conn, "gnomad") else None
    if is_populated(conn, "clinvar"):
        sources.append(ClinVar(conn, frequencies=frequencies))
    if is_populated(conn, "alphamissense"):
        sources.append(AlphaMissense(conn))
    if is_populated(conn, "gwas"):
        sources.append(GwasCatalog(conn))
    if is_populated(conn, "cpic"):
        sources.append(Cpic(conn))
    if snpedia_report:
        sources.append(SNPedia(snpedia_report))
    return sources


def _load_panel(path: str) -> set[str]:
    """Read an rsID panel: one per line, '#' comments and blanks ignored."""
    wanted = set()
    for line in Path(path).expanduser().read_text(encoding="utf-8").splitlines():
        token = line.split("#", 1)[0].strip()
        if token:
            wanted.add(token.lower())
    return wanted


def cmd_report(args: argparse.Namespace) -> int:
    path = Path(args.input).expanduser()
    sample = parse(path)
    print(
        f"Read {len(sample):,} variants from {sample.source_format} "
        f"(build {sample.build or 'unknown'})",
        file=sys.stderr,
    )
    for warning in sample.warnings:
        print(f"  warning: {warning}", file=sys.stderr)

    conn = connect(args.data_dir)

    # Check the declared build against the coordinates before anything joins on
    # a position. Requires both ClinVar builds indexed; harmless without them.
    if not args.no_build_check:
        before = sample.build
        sample = verify_build(sample, conn)
        if sample.build != before:
            print(f"  build resolved to GRCh{sample.build}", file=sys.stderr)

    # A Promethease report annotates itself; otherwise SNPedia needs one supplied.
    snpedia_report = args.snpedia or (str(path) if sample.source_format == "promethease" else None)
    sources = _build_sources(conn, snpedia_report)

    if not sources:
        print(
            "No annotation sources available. Run 'allele db update' first, "
            "or pass --snpedia with a Promethease report.",
            file=sys.stderr,
        )
        return 1
    print(f"Sources: {', '.join(s.name for s in sources)}", file=sys.stderr)

    if args.panel:
        wanted = _load_panel(args.panel)
        before = len(sample.calls)
        sample.calls = [
            c for c in sample.calls if c.rsid.lower() in wanted
        ]
        print(
            f"Panel filter: {len(sample.calls):,} of {before:,} variants kept "
            f"({len(wanted):,} requested)",
            file=sys.stderr,
        )

    report = annotate(sample, sources, materiality=not args.all)
    summary = report.summary()
    print(
        f"{len(report.credible):,} findings "
        f"({summary['clinical']:,} clinical, {summary['traits']:,} traits), "
        f"{summary['likely_artifacts']:,} demoted as likely artifacts",
        file=sys.stderr,
    )

    suffix = {"html": ".report.html", "json": ".report.json", "vcf": ".report.vcf"}
    output = Path(args.output) if args.output else path.with_suffix(suffix[args.format])
    if args.format == "json":
        body = to_json(report)
    elif args.format == "vcf":
        body = to_vcf(report)
    else:
        body = render(report, title=args.title)
    output.write_text(body, encoding="utf-8")
    print(f"Wrote {output}", file=sys.stderr)
    print(output)
    return 0


def cmd_diff(args: argparse.Namespace) -> int:
    import json

    before = json.loads(Path(args.before).expanduser().read_text(encoding="utf-8"))
    after = json.loads(Path(args.after).expanduser().read_text(encoding="utf-8"))
    result = diff(before, after)
    counts = result["counts"]
    print(
        f"{counts['added']} added, {counts['removed']} removed, "
        f"{counts['changed']} reclassified since {result['before']}"
    )
    for item in result["changed"][:40]:
        genes = ", ".join(item["genes"][:2]) if item["genes"] else ""
        print(f"  {item['rsid']:14s} {genes:20s} {item['changes']}")
    for finding in result["added"][:20]:
        print(f"  + {finding['rsid']:14s} {', '.join(finding['genes'][:2])}")
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    sample = parse(Path(args.input).expanduser())
    print(f"format   {sample.source_format}")
    print(f"build    {sample.build or 'not declared'}")
    for key, value in sample.stats().items():
        print(f"{key:9s} {value:,}")
    for warning in sample.warnings:
        print(f"warning  {warning}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="allele",
        description="Turn a consumer DNA export into an annotated report, offline.",
    )
    parser.add_argument("--version", action="version", version=f"allele {__version__}")
    parser.add_argument("--data-dir", help="where databases live (default: ~/.local/share/allele)")
    subparsers = parser.add_subparsers(dest="command", required=True)

    db = subparsers.add_parser("db", help="manage annotation databases")
    db.add_argument("action", choices=["update", "status", "clean", "path"])
    db.add_argument("--build", type=int, default=37, choices=[37, 38])
    db.add_argument("--only", choices=["clinvar", "gwas", "cpic", "gnomad", "alphamissense"])
    db.set_defaults(func=cmd_db)

    report = subparsers.add_parser("report", help="annotate a DNA export and write HTML")
    report.add_argument("input")
    report.add_argument("-o", "--output")
    report.add_argument("--snpedia", help="a Promethease report to use for SNPedia curation")
    report.add_argument("--title", default="Allele report")
    report.add_argument("--all", action="store_true", help="skip the materiality filter")
    report.add_argument("-f", "--format", choices=["html", "json", "vcf"], default="html")
    report.add_argument("--panel", help="file of rsIDs to restrict the report to")
    report.add_argument("--no-build-check", action="store_true",
                        help="trust the header build instead of verifying against coordinates")
    report.set_defaults(func=cmd_report)

    diff_cmd = subparsers.add_parser("diff", help="compare two JSON reports")
    diff_cmd.add_argument("before")
    diff_cmd.add_argument("after")
    diff_cmd.set_defaults(func=cmd_diff)

    info = subparsers.add_parser("info", help="describe a file without annotating it")
    info.add_argument("input")
    info.set_defaults(func=cmd_info)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
