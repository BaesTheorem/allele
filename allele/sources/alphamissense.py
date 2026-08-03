"""AlphaMissense pathogenicity predictions for missense variants.

ClinVar only knows what submitters have reviewed. Most missense variants have
never been classified by anyone, and for those a computational prediction is
the only signal available. AlphaMissense scores essentially all of them.

A prediction is not a classification, and the report says so wherever one
appears. It is evidence about a variant nobody has assessed, not a verdict.

Licensing: CC BY-NC-SA 4.0, the same non-commercial, share-alike terms as
SNPedia. It is therefore downloaded by the user and never redistributed, and
never enters the browser bundle. GRCh38 only, matching the file DeepMind
publishes.
"""

from __future__ import annotations

import gzip
import sqlite3
import urllib.request
from pathlib import Path
from typing import Iterator

from ..model import Call
from ..store import record_provenance, replace_table
from .base import Annotation, SourceInfo

URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
LICENSE = "AlphaMissense (DeepMind), CC BY-NC-SA 4.0. Not redistributed."
BUILD = 38

# DeepMind's own thresholds, from the AlphaMissense paper's class labels.
LIKELY_PATHOGENIC = 0.564
LIKELY_BENIGN = 0.34


def download(directory: Path, progress=None) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "AlphaMissense_hg38.tsv.gz"
    if target.is_file():
        return target
    if progress:
        progress("downloading AlphaMissense (0.6 GB)")
    urllib.request.urlretrieve(URL, target)
    return target


def _iter_records(path: Path, wanted: dict[str, set[int]]) -> Iterator[tuple]:
    """Stream the file, keeping only positions we might ever look up.

    The full table is ~71M rows. Filtering to positions ClinVar flags as
    notable keeps the index small enough to query instantly, at the cost of
    having nothing to say about variants outside that set.
    """
    seen: set[tuple[str, int, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0][3:] if fields[0].startswith("chr") else fields[0]
            positions = wanted.get(chrom)
            if not positions:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            if pos not in positions:
                continue
            ref, alt = fields[2], fields[3]
            key = (chrom, pos, alt)
            if key in seen:
                continue
            seen.add(key)
            try:
                score = float(fields[8])
            except ValueError:
                continue
            yield (chrom, pos, BUILD, ref, alt, score, fields[9])


def build_index(conn: sqlite3.Connection, directory: Path, progress=None) -> int:
    from .clinvar import is_notable

    wanted: dict[str, set[int]] = {}
    query = "SELECT chrom, pos, significance FROM clinvar WHERE build = ? AND pos IS NOT NULL"
    for row in conn.execute(query, (BUILD,)):
        if is_notable(row["significance"]):
            wanted.setdefault(row["chrom"], set()).add(row["pos"])
    if not wanted:
        raise RuntimeError("no GRCh38 ClinVar positions indexed; run 'allele db update --build 38'")

    path = download(directory, progress)
    if progress:
        progress(f"indexing AlphaMissense over {sum(len(v) for v in wanted.values()):,} positions")
    columns = ["chrom", "pos", "build", "ref", "alt", "score", "class"]
    count = replace_table(conn, "alphamissense", columns, _iter_records(path, wanted))
    record_provenance(conn, "alphamissense", "hg38", URL, LICENSE, count)
    return count


class AlphaMissense:
    """Computational pathogenicity for variants nobody has classified."""

    name = "alphamissense"

    def __init__(self, conn: sqlite3.Connection):
        self.conn = conn

    def info(self) -> SourceInfo:
        row = self.conn.execute(
            "SELECT * FROM provenance WHERE source='alphamissense'"
        ).fetchone()
        if not row:
            return SourceInfo(name="alphamissense", license=LICENSE)
        return SourceInfo(
            name="alphamissense",
            release=row["release"],
            downloaded=row["downloaded"],
            url=row["url"],
            license=row["license"],
            record_count=row["record_count"],
            notes=[
                "prediction, not a classification",
                "GRCh38 missense variants only",
            ],
        )

    def score(self, chrom: str | None, pos: int | None, alt: str) -> tuple[float, str] | None:
        if not chrom or not pos:
            return None
        row = self.conn.execute(
            "SELECT score, class FROM alphamissense WHERE chrom=? AND pos=? AND alt=? LIMIT 1",
            (chrom, pos, alt),
        ).fetchone()
        if not row:
            return None
        return row["score"], row["class"]

    def lookup(self, call: Call) -> list[Annotation]:
        """Only speaks where ClinVar is silent.

        Where an expert panel has already classified a variant, a prediction
        adds nothing and would only muddy the picture. This is here for the
        long tail nobody has reviewed.
        """
        if not call.chrom or not call.pos or not call.usable:
            return []
        rows = self.conn.execute(
            "SELECT * FROM alphamissense WHERE chrom=? AND pos=?", (call.chrom, call.pos)
        ).fetchall()

        out: list[Annotation] = []
        for row in rows:
            alt = row["alt"]
            if alt not in call.genotype:
                continue
            if row["score"] < LIKELY_PATHOGENIC:
                continue  # only surface the predicted-pathogenic tail
            zygosity = (
                "homozygous" if call.genotype.count(alt) == 2
                else "hemizygous" if len(call.genotype) == 1
                else "heterozygous"
            )
            out.append(
                Annotation(
                    rsid=call.rsid,
                    source="alphamissense",
                    category="predicted",
                    title=f"predicted {row['class'].replace('_', ' ')} ({row['score']:.2f})",
                    significance=row["class"].replace("_", " "),
                    zygosity=zygosity,
                    genotype=call.genotype,
                    flags=(
                        "computational prediction, not a classification by anyone. "
                        "AlphaMissense scores structure and conservation; it has no "
                        "clinical evidence behind it.",
                    ),
                )
            )
        return out
