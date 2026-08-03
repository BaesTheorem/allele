"""Read consumer DNA exports into normalized `Call` records.

Formats are detected from file content, never from the extension, because
vendors are careless about it (Living DNA ships tab-delimited data in a file
named .csv). Everything streams: these files run to hundreds of megabytes
uncompressed and there is no reason to hold one in memory.

Strand convention: every consumer chip export listed here reports genotypes on
the plus strand of the stated reference assembly, which is also what VCF
REF/ALT use. So raw exports need no flipping to match ClinVar. Promethease
reports are the exception and carry explicit per-SNP orientation, which
`promethease.py` preserves.
"""

from __future__ import annotations

import csv
import gzip
import io
import re
import zipfile
from pathlib import Path
from typing import Iterable, Iterator

from .model import (
    Call,
    Sample,
    normalize_chromosome,
    normalize_genotype,
    normalize_ploidy,
)

# "build 37", "build 36.3", "GRCh38", "human assembly build 37"
_BUILD_RE = re.compile(r"(?:build|grch)\s*([0-9]{2})", re.IGNORECASE)
_HUMAN_BUILD = {"36": 36, "37": 37, "38": 38, "19": 37, "18": 36}


def _sniff_build(text: str) -> int | None:
    """Pull a genome build out of a comment header."""
    for match in _BUILD_RE.finditer(text):
        build = _HUMAN_BUILD.get(match.group(1))
        if build:
            return build
    return None


def open_text(path: Path) -> tuple[io.TextIOBase, str]:
    """Open a report whether it is plain, gzipped, or a single-member zip."""
    with open(path, "rb") as probe:
        magic = probe.read(4)

    if magic[:2] == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace"), "gzip"

    if magic[:2] == b"PK":
        archive = zipfile.ZipFile(path)
        members = [n for n in archive.namelist() if not n.endswith("/")]
        if not members:
            raise ValueError(f"{path.name}: zip archive is empty")
        # Vendor zips contain one data file plus occasional metadata.
        members.sort(key=lambda n: archive.getinfo(n).file_size, reverse=True)
        stream = archive.open(members[0])
        return io.TextIOWrapper(stream, encoding="utf-8", errors="replace"), f"zip:{members[0]}"

    return open(path, "r", encoding="utf-8", errors="replace"), "text"


def _read_header(handle: io.TextIOBase, max_lines: int = 60) -> tuple[list[str], list[str]]:
    """Consume leading comment lines. Returns (comments, pending data lines)."""
    comments: list[str] = []
    pending: list[str] = []
    for _ in range(max_lines):
        line = handle.readline()
        if not line:
            break
        if line.startswith("#"):
            comments.append(line.rstrip("\n"))
            continue
        pending.append(line)
        break
    return comments, pending


def detect_format(comments: list[str], first_data_line: str) -> str:
    """Identify the vendor from the header text and the shape of row one."""
    header = " ".join(comments).lower()
    first = first_data_line.strip()
    lower_first = first.lower()

    if "23andme" in header:
        return "23andme"
    if "ancestrydna" in header or "ancestry.com" in header:
        return "ancestrydna"
    if "myheritage" in header:
        return "myheritage"
    if "living dna" in header or "livingdna" in header:
        return "livingdna"
    # Most specific first: FamFinder headers also contain "family tree dna".
    if "famfinder" in header:
        return "ftdna-famfinder"
    if "myhappygenes" in header or "tempus" in header:
        return "myhappygenes"
    if "ftdna" in header or "family tree dna" in header or "familytreedna" in header:
        return "ftdna"

    # No usable header. Fall back to column shape.
    if lower_first.startswith("rsid,chromosome,position,result"):
        return "ftdna"  # also MyHeritage; both parse identically
    if lower_first.startswith("rsid\tchromosome\tposition\tallele1"):
        # FamFinder uses the same header; the vendor header line disambiguates.
        return "ftdna-famfinder" if "famfinder" in header else "ancestrydna"
    if lower_first.startswith("rsid\tchromosome\tposition\tresult"):
        return "ftdna-illumina"
    if lower_first.startswith("rsid"):
        return "generic"

    fields = re.split(r"[\t,]", first)
    if len(fields) == 4 and fields[0].lower().startswith(("rs", "i")):
        return "generic"
    if len(fields) == 5 and fields[0].lower().startswith(("rs", "i")):
        return "ancestrydna"
    return "unknown"


def _split(line: str) -> list[str]:
    """Split on tab, or comma if the line has no tabs."""
    stripped = line.rstrip("\n\r")
    parts = stripped.split("\t") if "\t" in stripped else next(csv.reader([stripped]))
    return [p.strip().strip('"') for p in parts]


def _iter_chip_rows(
    lines: Iterable[str], fmt: str, build: int | None, warnings: list[str]
) -> Iterator[Call]:
    """Parse the flat 4- or 5-column chip layouts shared by every vendor."""
    seen_header = False
    malformed = 0

    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        fields = _split(line)
        if len(fields) < 4:
            malformed += 1
            continue

        # Some exports repeat a column header as the first data row.
        if not seen_header and fields[0].lower() in ("rsid", "rs id", "snp"):
            seen_header = True
            continue

        rsid = fields[0]
        chrom = normalize_chromosome(fields[1])
        raw_pos = fields[2]

        # 4 columns: concatenated genotype. 5 columns: separate alleles.
        raw_genotype = fields[3] if len(fields) == 4 else fields[3] + fields[4]

        try:
            pos = int(raw_pos)
        except (TypeError, ValueError):
            pos = None

        yield Call(
            rsid=rsid,
            chrom=chrom,
            pos=pos,
            genotype=normalize_ploidy(chrom, normalize_genotype(raw_genotype)),
            build=build,
        )

    if malformed:
        warnings.append(f"{malformed} malformed line(s) skipped")


def _iter_vcf_rows(handle: io.TextIOBase, warnings: list[str]) -> tuple[list[Call], int | None]:
    """Parse a single-sample VCF into calls, resolving GT indices to alleles."""
    calls: list[Call] = []
    build: int | None = None
    skipped_multisample = False
    gvcf_blocks = 0

    for line in handle:
        if line.startswith("##"):
            if build is None and ("reference=" in line or "assembly=" in line):
                build = _sniff_build(line)
            continue
        if line.startswith("#CHROM"):
            columns = line.rstrip("\n").split("\t")
            if len(columns) > 10 and not skipped_multisample:
                warnings.append(
                    f"multi-sample VCF: using the first sample ({columns[9]}) of {len(columns) - 9}"
                )
                skipped_multisample = True
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue

        chrom = normalize_chromosome(fields[0])
        rsid = fields[2]
        ref, alt_field = fields[3], fields[4]

        # gVCF reference blocks (<NON_REF> with an END= span) assert only that a
        # region matched the reference. They are not variant calls.
        if "<NON_REF>" in alt_field or "<*>" in alt_field:
            gvcf_blocks += 1
            continue
        fmt_keys = fields[8].split(":")
        sample_values = fields[9].split(":")

        if "GT" not in fmt_keys:
            continue
        gt = sample_values[fmt_keys.index("GT")]
        if gt in (".", "./.", ".|."):
            continue

        alleles = [ref] + alt_field.split(",")
        bases = []
        ok = True
        for index in re.split(r"[/|]", gt):
            if index == ".":
                ok = False
                break
            try:
                allele = alleles[int(index)]
            except (ValueError, IndexError):
                ok = False
                break
            # Only single-base alleles map onto the chip-style genotype model.
            if len(allele) != 1:
                ok = False
                break
            bases.append(allele)
        if not ok or not bases:
            continue

        try:
            pos = int(fields[1])
        except ValueError:
            pos = None

        calls.append(
            Call(
                rsid=rsid if rsid and rsid != "." else "",
                chrom=chrom,
                pos=pos,
                genotype=normalize_ploidy(chrom, normalize_genotype("".join(bases))),
                build=build,
            )
        )

    if gvcf_blocks:
        warnings.append(
            f"{gvcf_blocks:,} gVCF reference block(s) skipped; they assert reference "
            "match over a span rather than a genotype call"
        )
    return calls, build


def parse(path: str | Path) -> Sample:
    """Read any supported export into a `Sample`.

    Promethease HTML reports are handled by `promethease.parse_report`; this
    function dispatches to it so callers have one entry point.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)

    handle, container = open_text(path)
    try:
        head = handle.read(4096)
        if "<html" in head[:2048].lower() or "promethease" in head.lower():
            from .promethease import parse_report

            handle.close()
            return parse_report(path)

        handle.seek(0) if container == "text" else None
        if container != "text":
            # TextIOWrapper over a compressed stream is not seekable; reopen.
            handle.close()
            handle, container = open_text(path)

        if head.lstrip().startswith("##fileformat=VCF"):
            calls, build = _iter_vcf_rows(handle, warnings := [])
            sample = Sample(
                calls=calls,
                source_format="vcf",
                build=build,
                path=str(path),
                warnings=warnings,
            )
            if build is None:
                sample.warnings.append(
                    "VCF does not declare a reference assembly; positions cannot be verified"
                )
            return sample

        comments, pending = _read_header(handle)
        first_data_line = pending[0] if pending else ""
        fmt = detect_format(comments, first_data_line)
        build = _sniff_build(" ".join(comments))
        warnings: list[str] = []

        if fmt == "unknown":
            warnings.append(
                "unrecognized layout; parsed as a generic rsid/chromosome/position/genotype table"
            )
            fmt = "generic"

        lines = iter(pending + handle.readlines())
        calls = list(_iter_chip_rows(lines, fmt, build, warnings))

        if build is None:
            # Every vendor here has shipped build 37 for years, but guessing
            # silently is how positions end up quietly wrong. Say so instead.
            warnings.append(
                "no genome build declared in the header; "
                "position-based annotation will be skipped (rsID matching still works)"
            )

        return Sample(
            calls=calls,
            source_format=fmt,
            build=build,
            path=str(path),
            warnings=warnings,
        )
    finally:
        if not handle.closed:
            handle.close()
