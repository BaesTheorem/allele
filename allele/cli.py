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
from .parsers import parse
from .report import render
from .store import DEFAULT_DIR, connect, finalize, is_populated, provenance


def _progress(message: str) -> None:
    print(f"  {message}", file=sys.stderr, flush=True)


def cmd_db(args: argparse.Namespace) -> int:
    from .sources import clinvar, gwas

    directory = Path(args.data_dir) if args.data_dir else DEFAULT_DIR
    downloads = directory / "downloads"
    conn = connect(directory)

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
    finalize(conn)
    print("Done.", file=sys.stderr)
    return 0


def _build_sources(conn, snpedia_report: str | None) -> list:
    from .sources import ClinVar, GwasCatalog, SNPedia

    sources: list = []
    if is_populated(conn, "clinvar"):
        sources.append(ClinVar(conn))
    if is_populated(conn, "gwas"):
        sources.append(GwasCatalog(conn))
    if snpedia_report:
        sources.append(SNPedia(snpedia_report))
    return sources


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

    report = annotate(sample, sources, materiality=not args.all)
    summary = report.summary()
    print(
        f"{len(report.credible):,} findings "
        f"({summary['clinical']:,} clinical, {summary['traits']:,} traits), "
        f"{summary['likely_artifacts']:,} demoted as likely artifacts",
        file=sys.stderr,
    )

    output = Path(args.output) if args.output else path.with_suffix(".report.html")
    output.write_text(render(report, title=args.title), encoding="utf-8")
    print(f"Wrote {output}", file=sys.stderr)
    print(output)
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
    db.add_argument("action", choices=["update", "status"])
    db.add_argument("--build", type=int, default=37, choices=[37, 38])
    db.add_argument("--only", choices=["clinvar", "gwas"])
    db.set_defaults(func=cmd_db)

    report = subparsers.add_parser("report", help="annotate a DNA export and write HTML")
    report.add_argument("input")
    report.add_argument("-o", "--output")
    report.add_argument("--snpedia", help="a Promethease report to use for SNPedia curation")
    report.add_argument("--title", default="Allele report")
    report.add_argument("--all", action="store_true", help="skip the materiality filter")
    report.set_defaults(func=cmd_report)

    info = subparsers.add_parser("info", help="describe a file without annotating it")
    info.add_argument("input")
    info.set_defaults(func=cmd_info)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
