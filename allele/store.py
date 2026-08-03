"""Local SQLite store for downloaded annotation databases.

One file, one table per source, plus a provenance table so a report can state
which release each finding came from and how old it is. Everything is written
in a single transaction per source and indexed on rsID, which is the only key
a consumer chip gives us.
"""

from __future__ import annotations

import sqlite3
from datetime import date
from pathlib import Path
from typing import Iterable, Iterator

DEFAULT_DIR = Path.home() / ".local" / "share" / "allele"

SCHEMA = """
CREATE TABLE IF NOT EXISTS provenance (
    source      TEXT PRIMARY KEY,
    release     TEXT,
    downloaded  TEXT,
    url         TEXT,
    license     TEXT,
    record_count INTEGER
);

CREATE TABLE IF NOT EXISTS clinvar (
    rsid         TEXT NOT NULL,
    chrom        TEXT,
    pos          INTEGER,
    build        INTEGER NOT NULL,
    ref          TEXT,
    alt          TEXT,
    significance TEXT,
    review       TEXT,
    conditions   TEXT,
    genes        TEXT,
    frequency    REAL
);

CREATE TABLE IF NOT EXISTS gnomad (
    chrom     TEXT NOT NULL,
    pos       INTEGER NOT NULL,
    build     INTEGER NOT NULL,
    ref       TEXT,
    alt       TEXT,
    af        REAL,
    af_grpmax REAL
);

CREATE TABLE IF NOT EXISTS alphamissense (
    chrom     TEXT NOT NULL,
    pos       INTEGER NOT NULL,
    build     INTEGER NOT NULL,
    ref       TEXT,
    alt       TEXT,
    score     REAL,
    class     TEXT
);

CREATE TABLE IF NOT EXISTS cpic (
    rsid         TEXT NOT NULL,
    gene         TEXT,
    variant_name TEXT,
    drug         TEXT,
    level        TEXT,
    guideline    TEXT,
    url          TEXT
);

CREATE TABLE IF NOT EXISTS gwas (
    rsid        TEXT NOT NULL,
    risk_allele TEXT,
    trait       TEXT,
    mapped_trait TEXT,
    p_value     REAL,
    effect_size TEXT,
    genes       TEXT,
    pubmed_id   TEXT,
    first_author TEXT,
    year        TEXT,
    journal     TEXT,
    study_title TEXT,
    accession   TEXT,
    cohort      TEXT,
    cohort_size INTEGER,
    ancestries  TEXT,
    risk_freq   REAL
);
"""

INDEXES = """
CREATE INDEX IF NOT EXISTS idx_clinvar_rsid ON clinvar(rsid);
CREATE INDEX IF NOT EXISTS idx_clinvar_rsid_build ON clinvar(rsid, build);
CREATE INDEX IF NOT EXISTS idx_gwas_rsid ON gwas(rsid);
CREATE INDEX IF NOT EXISTS idx_cpic_rsid ON cpic(rsid);
CREATE INDEX IF NOT EXISTS idx_gnomad_pos ON gnomad(chrom, pos, alt);
CREATE INDEX IF NOT EXISTS idx_am_pos ON alphamissense(chrom, pos, alt);
"""


def connect(directory: Path | str | None = None) -> sqlite3.Connection:
    """Open (creating if needed) the local database."""
    base = Path(directory) if directory else DEFAULT_DIR
    base.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(base / "annotations.sqlite")
    conn.row_factory = sqlite3.Row
    conn.executescript(SCHEMA)
    return conn


def finalize(conn: sqlite3.Connection) -> None:
    """Build indexes after bulk load. Much faster than indexing as we go."""
    conn.executescript(INDEXES)
    conn.commit()


def replace_table(
    conn: sqlite3.Connection,
    table: str,
    columns: list[str],
    rows: Iterable[tuple],
    where: tuple[str, tuple] | None = None,
) -> int:
    """Swap a table's contents for a fresh load, counting what went in.

    `where` scopes the delete, so ClinVar can hold GRCh37 and GRCh38 side by
    side and rebuilding one build does not wipe the other.
    """
    placeholders = ",".join("?" * len(columns))
    if where:
        conn.execute(f"DELETE FROM {table} WHERE {where[0]}", where[1])
    else:
        conn.execute(f"DELETE FROM {table}")
    count = 0

    def counted() -> Iterator[tuple]:
        nonlocal count
        for row in rows:
            count += 1
            yield row

    conn.executemany(
        f"INSERT INTO {table} ({','.join(columns)}) VALUES ({placeholders})", counted()
    )
    conn.commit()
    return count


def record_provenance(
    conn: sqlite3.Connection,
    source: str,
    release: str | None,
    url: str,
    license_text: str,
    record_count: int,
) -> None:
    conn.execute(
        "INSERT OR REPLACE INTO provenance "
        "(source, release, downloaded, url, license, record_count) VALUES (?,?,?,?,?,?)",
        (source, release, date.today().isoformat(), url, license_text, record_count),
    )
    conn.commit()


def provenance(conn: sqlite3.Connection) -> dict[str, sqlite3.Row]:
    return {row["source"]: row for row in conn.execute("SELECT * FROM provenance")}


def is_populated(conn: sqlite3.Connection, table: str) -> bool:
    try:
        return conn.execute(f"SELECT 1 FROM {table} LIMIT 1").fetchone() is not None
    except sqlite3.OperationalError:
        return False
