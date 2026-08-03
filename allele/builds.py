"""Detect a file's genome build from its coordinates, not just its header.

Headers lie. Real exports have been observed declaring "build 37" while
carrying GRCh38 coordinates, and a header taken on trust silently shifts every
position by however much the assemblies differ.

The check is cheap because we already hold both ClinVar builds: any rsID whose
GRCh37 and GRCh38 positions differ is a discriminating probe. Look up a sample
of the file's own variants and see which assembly its positions agree with.
"""

from __future__ import annotations

import sqlite3

from .model import Sample

# How many probes to check, and how decisive the winner must be. A file with
# no clear majority gets no verdict rather than a coin flip.
SAMPLE_SIZE = 400
MIN_PROBES = 20
WIN_RATIO = 0.8


def load_probes(conn: sqlite3.Connection, limit: int = 20000) -> dict[str, dict[int, int]]:
    """rsIDs whose position differs between builds, with both positions."""
    query = """
        SELECT a.rsid, a.pos AS pos37, b.pos AS pos38
        FROM clinvar a
        JOIN clinvar b ON a.rsid = b.rsid AND a.build = 37 AND b.build = 38
        WHERE a.pos IS NOT NULL AND b.pos IS NOT NULL AND a.pos != b.pos
        LIMIT ?
    """
    probes: dict[str, dict[int, int]] = {}
    for row in conn.execute(query, (limit,)):
        probes.setdefault(row["rsid"], {37: row["pos37"], 38: row["pos38"]})
    return probes


def detect(sample: Sample, conn: sqlite3.Connection) -> tuple[int | None, dict]:
    """Infer the build from coordinates. Returns (build, evidence)."""
    probes = load_probes(conn)
    if not probes:
        return None, {"reason": "no dual-build reference data indexed"}

    votes = {37: 0, 38: 0}
    checked = 0
    for call in sample.calls:
        if checked >= SAMPLE_SIZE:
            break
        if not call.pos or not call.rsid:
            continue
        expected = probes.get(call.rsid)
        if not expected:
            continue
        checked += 1
        for build, position in expected.items():
            if call.pos == position:
                votes[build] += 1

    total = votes[37] + votes[38]
    evidence = {"probes_checked": checked, "matched": total, "votes": dict(votes)}

    if total < MIN_PROBES:
        evidence["reason"] = "too few positions matched a dual-build reference variant"
        return None, evidence

    winner = max(votes, key=votes.get)
    if votes[winner] / total < WIN_RATIO:
        evidence["reason"] = "no assembly won decisively"
        return None, evidence

    evidence["confidence"] = round(votes[winner] / total, 3)
    return winner, evidence


def verify(sample: Sample, conn: sqlite3.Connection) -> Sample:
    """Check the declared build against the coordinates and reconcile.

    A header that disagrees with the data is not a tie: the coordinates are the
    thing annotation actually joins on, so they win, and the disagreement is
    recorded rather than smoothed over.
    """
    detected, evidence = detect(sample, conn)
    if detected is None:
        if sample.build is None:
            sample.warnings.append(
                "Genome build could not be determined from the header or the "
                f"coordinates ({evidence.get('reason', 'inconclusive')}). "
                "Annotation falls back to rsID matching."
            )
        return sample

    if sample.build is None:
        sample.build = detected
        sample.calls = [
            type(c)(c.rsid, c.chrom, c.pos, c.genotype, detected) for c in sample.calls
        ]
        sample.warnings.append(
            f"No build declared in the header; inferred GRCh{detected} from coordinates "
            f"({evidence['votes'][detected]}/{evidence['matched']} probes agree)."
        )
    elif sample.build != detected:
        sample.warnings.append(
            f"Header declares build {sample.build} but the coordinates match "
            f"GRCh{detected} ({evidence['votes'][detected]}/{evidence['matched']} probes). "
            "Using the coordinates, since that is what annotation joins on."
        )
        sample.build = detected
        sample.calls = [
            type(c)(c.rsid, c.chrom, c.pos, c.genotype, detected) for c in sample.calls
        ]
    return sample
