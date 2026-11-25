"""Check alignments against the genomes they were made from."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator
from typing import NamedTuple

from .blocks import GAP, Row, reverse_complement
from .fasta import Genomes
from .mafio import MafError, MafReader

HEADER = ("line", "src", "start", "size", "strand", "problem", "detail")

# Shorter rows match by chance, so no explanation is attempted for them
MIN_EXPLAINED = 12
MAX_SHIFT = 10


class Problem(NamedTuple):
    line: int | str
    src: str
    start: int | str
    size: int | str
    strand: str
    problem: str
    detail: str


def _fetch(genomes: Genomes, row: Row, start: int, strand: str) -> str:
    """Bases of the row's length from its source, at a start on a given strand."""
    fasta, name = genomes.find(row.src)
    if strand == "+":
        return fasta.fetch(name, start, start + row.size)
    length = fasta.records[name].length
    return reverse_complement(fasta.fetch(name, length - start - row.size, length - start))


def _differences(bases: str, truth: str) -> int:
    """Number of differing bases. N in the alignment stands for any base."""
    if bases == truth:
        return 0
    if len(bases) != len(truth):
        return max(len(bases), len(truth))
    return sum(a != b and a != "N" for a, b in zip(bases, truth))


def _explain(genomes: Genomes, row: Row, bases: str, length: int) -> tuple[str, str] | None:
    """A simple mistake that accounts for the bases of a row, if there is one."""
    if row.size < MIN_EXPLAINED:
        return None
    other = "-" if row.strand == "+" else "+"
    mirrored = length - row.end

    def fits(start: int, strand: str) -> bool:
        return 0 <= start <= length - row.size and not _differences(
            bases, _fetch(genomes, row, start, strand)
        )

    if fits(row.start, other):
        return "strand", f"bases match strand {other} at the same start"
    if fits(mirrored, row.strand):
        return "start", f"start counts from the other end of the sequence, expected {mirrored}"
    if fits(mirrored, other):
        return "orientation", "bases are the reverse complement"
    for shift in sorted(range(-MAX_SHIFT, MAX_SHIFT + 1), key=abs):
        if fits(row.start + shift, row.strand):
            return "start", f"shifted by {-shift:+d}, expected {row.start + shift}"
    return None


def check_row(row: Row, genomes: Genomes) -> tuple[str, str] | None:
    """The problem of a row and its details, or None. The source must be known."""
    fasta, name = genomes.find(row.src)
    length = fasta.records[name].length
    if row.src_size != length:
        return "srcSize", f"the sequence has {length} bases, not {row.src_size}"
    if row.start < 0 or row.end > length:
        return "range", f"outside of the sequence of {length} bases"
    bases = row.text.replace(GAP, "").upper().replace("U", "T")
    differences = _differences(bases, _fetch(genomes, row, row.start, row.strand))
    if not differences:
        return None
    return _explain(genomes, row, bases, length) or (
        "bases", f"{differences} of {row.size} differ"
    )


def verify(
    reader: MafReader, genomes: Genomes, stats: Counter, unknown: Counter
) -> Iterator[Problem]:
    """Check every row and yield the problems.

    Malformed lines are reported as problems instead of ending the check. Rows
    without a sequence are counted in ``unknown``, by source.
    """
    malformed: list[MafError] = []
    reader.on_error = malformed.append

    def report_malformed() -> Iterator[Problem]:
        for error in malformed:
            stats["checked"] += 1
            stats["problems"] += 1
            yield Problem(error.line, ".", ".", ".", ".", "format", error.reason)
        malformed.clear()

    for block in reader:
        yield from report_malformed()
        for row in block.rows:
            if genomes.find(row.src) is None:
                unknown[row.src] += 1
                continue
            stats["checked"] += 1
            found = check_row(row, genomes)
            if found is not None:
                stats["problems"] += 1
                yield Problem(reader.block_line, row.src, row.start, row.size, row.strand, *found)
    yield from report_malformed()
