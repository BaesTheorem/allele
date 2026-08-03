"""gnomAD population allele frequencies, fetched by targeted range request.

Frequency is the single most useful number for deciding whether a consumer
array's "pathogenic" call is real. Without it the tool has to infer rarity from
a variant's absence in older, smaller cohorts, which is defensible but weak.

The full gnomAD release is roughly 100 GB, which is not a reasonable download.
It is also tabix-indexed and served from S3 with range-request support, so only
the blocks covering variants we actually care about need fetching. Merging the
notable ClinVar positions into windows brings that to a few thousand requests.

GRCh38 only, which is what gnomAD v4 publishes. Samples on GRCh37 get no
frequency from here, and the cohort-absence rule still applies to them.
"""

from __future__ import annotations

import sqlite3
from concurrent.futures import ThreadPoolExecutor
from typing import Iterator

from ..model import Call
from ..store import record_provenance, replace_table
from .base import Annotation, SourceInfo

BASE = (
    "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/exomes/"
    "gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"
)
LICENSE = "gnomAD, ODbL v1.0 (Broad Institute)"
BUILD = 38

# Merging nearby target positions into one window trades a little wasted
# parsing for far fewer round trips; 200 kb keeps requests in the low thousands.
# Large windows pull enough data that S3 connections drop mid-stream and
# htslib reports a truncated file, so keep each request modest and retry.
WINDOW_MERGE_GAP = 20_000
MAX_RETRIES = 3
# Purely I/O bound: every worker is waiting on S3, so threads are the right
# tool, and the chromosome is the natural unit of work.
WORKERS = 8
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]


def _windows(positions: list[int], gap: int = WINDOW_MERGE_GAP) -> Iterator[tuple[int, int]]:
    start = end = None
    for pos in positions:
        if start is None:
            start, end = pos, pos
        elif pos <= end + gap:
            end = pos
        else:
            yield start, end
            start, end = pos, pos
    if start is not None:
        yield start, end


def targets(conn: sqlite3.Connection) -> dict[str, list[int]]:
    """Positions worth fetching: notable ClinVar variants on GRCh38."""
    from .clinvar import is_notable

    by_chrom: dict[str, set[int]] = {}
    query = "SELECT chrom, pos, significance FROM clinvar WHERE build = ? AND pos IS NOT NULL"
    for row in conn.execute(query, (BUILD,)):
        if not is_notable(row["significance"]):
            continue
        by_chrom.setdefault(row["chrom"], set()).add(row["pos"])
    return {c: sorted(p) for c, p in by_chrom.items()}


def _extract_chromosome(chrom: str, positions: list[int], progress=None) -> list[tuple]:
    """Fetch every wanted position on one chromosome. Runs in its own thread."""
    import pysam

    wanted_set = set(positions)
    windows = list(_windows(positions))
    rows: list[tuple] = []
    try:
        handle = pysam.VariantFile(BASE.format(chrom=chrom))
    except Exception as exc:
        if progress:
            progress(f"  chr{chrom} unavailable: {str(exc)[:80]}")
        return rows

    found = failed = 0
    for number, (start, end) in enumerate(windows, 1):
        for attempt in range(MAX_RETRIES):
            try:
                for record in handle.fetch(f"chr{chrom}", max(0, start - 1), end + 1):
                    if record.pos not in wanted_set:
                        continue
                    # Only single-base substitutions map onto array genotypes.
                    if len(record.ref) != 1 or not record.alts:
                        continue
                    for index, alt in enumerate(record.alts):
                        if len(alt) != 1:
                            continue
                        frequency = _first(record.info.get("AF"), index)
                        grpmax = _first(record.info.get("AF_grpmax"), index)
                        if frequency is None and grpmax is None:
                            continue
                        rows.append((chrom, record.pos, BUILD, record.ref, alt, frequency, grpmax))
                        found += 1
                break
            except Exception:
                if attempt == MAX_RETRIES - 1:
                    failed += 1
                    break
                # A dropped stream leaves the handle unusable; reopen it.
                try:
                    handle.close()
                except Exception:
                    pass
                try:
                    handle = pysam.VariantFile(BASE.format(chrom=chrom))
                except Exception:
                    failed += 1
                    break
    try:
        handle.close()
    except Exception:
        pass
    if progress:
        note = f", {failed} window(s) failed" if failed else ""
        progress(f"  chr{chrom}: {found:,} frequencies from {len(windows)} windows{note}")
    return rows


def build_index(conn: sqlite3.Connection, progress=None, chromosomes=None, workers=WORKERS) -> int:
    """Fetch frequencies for the target positions and index them."""
    try:
        import pysam  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "gnomAD extraction needs pysam: pip install 'allele[gnomad]'"
        ) from exc

    wanted = targets(conn)
    if not wanted:
        raise RuntimeError("no GRCh38 ClinVar positions indexed; run 'allele db update --build 38'")

    todo = [(c, wanted[c]) for c in (chromosomes or CHROMOSOMES) if wanted.get(c)]
    if progress:
        total = sum(len(list(_windows(p))) for _, p in todo)
        progress(f"gnomAD: {total:,} windows across {len(todo)} chromosomes, {workers} workers")

    rows: list[tuple] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(lambda item: _extract_chromosome(item[0], item[1], progress), todo):
            rows.extend(result)

    columns = ["chrom", "pos", "build", "ref", "alt", "af", "af_grpmax"]
    count = replace_table(conn, "gnomad", columns, iter(rows))
    record_provenance(conn, "gnomad", "v4.1 exomes", BASE, LICENSE, count)
    return count


def _first(value, index: int):
    if value is None:
        return None
    if isinstance(value, (tuple, list)):
        if index < len(value):
            value = value[index]
        else:
            return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class Gnomad:
    """Frequency enrichment. Adds no findings of its own.

    This source exists to make other sources' findings judgeable: it supplies
    the population frequency the plausibility checks need. It never produces a
    finding by itself, because "this variant is common" is not a result.
    """

    name = "gnomad"

    def __init__(self, conn: sqlite3.Connection):
        self.conn = conn

    def info(self) -> SourceInfo:
        row = self.conn.execute("SELECT * FROM provenance WHERE source='gnomad'").fetchone()
        if not row:
            return SourceInfo(name="gnomad", license=LICENSE)
        return SourceInfo(
            name="gnomad",
            release=row["release"],
            downloaded=row["downloaded"],
            url=row["url"],
            license=row["license"],
            record_count=row["record_count"],
            notes=["GRCh38 only", "exomes; covers coding regions"],
        )

    def frequency(self, chrom: str | None, pos: int | None, alt: str) -> float | None:
        """Population frequency for one alternate allele, if gnomAD has it."""
        if not chrom or not pos:
            return None
        row = self.conn.execute(
            "SELECT af, af_grpmax FROM gnomad WHERE chrom=? AND pos=? AND alt=? LIMIT 1",
            (chrom, pos, alt),
        ).fetchone()
        if not row:
            return None
        return row["af"] if row["af"] is not None else row["af_grpmax"]

    def lookup(self, call: Call) -> list[Annotation]:
        return []
