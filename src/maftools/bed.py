"""BED input and interval lookup."""

from __future__ import annotations

import re
from bisect import bisect_left
from collections.abc import Iterable
from itertools import accumulate
from typing import NamedTuple

from .mafio import MafError

# e.g. NC_037638.1: the suffix is a version, not a chromosome name
_ACCESSION = re.compile(r"[A-Z]{1,6}_?\d{5,}\.\d+")


class Interval(NamedTuple):
    chrom: str
    start: int
    end: int
    name: str = ""
    strand: str = "."


def read_bed(
    handle: Iterable[str], name: str = "<bed>", stranded: bool = False
) -> list[Interval]:
    """Read chrom, start, end and, if present, name and strand.

    An invalid strand is an error if ``stranded``, and taken as unknown otherwise.
    """
    intervals = []
    for number, line in enumerate(handle, 1):
        # names may contain spaces if the columns are separated by tabs
        fields = [f.strip() for f in line.split("\t")] if "\t" in line else line.split()
        if not line.strip() or fields[0].startswith("#") or fields[0].split()[:1] in (["track"], ["browser"]):
            continue
        try:
            interval = Interval(fields[0], int(fields[1]), int(fields[2]), *fields[3:4])
        except (IndexError, ValueError):
            raise MafError(f"{name}:{number}: expected chrom, start and end") from None
        if not 0 <= interval.start <= interval.end:
            raise MafError(f"{name}:{number}: start must be >= 0 and <= end")
        if len(fields) == 5 and fields[4] in ("+", "-"):
            # not standard, but a score cannot be mistaken for a strand
            interval = interval._replace(strand=fields[4])
        elif len(fields) > 5:
            if fields[5] in ("+", "-"):
                interval = interval._replace(strand=fields[5])
            elif stranded and fields[5] != ".":
                raise MafError(f"{name}:{number}: invalid strand {fields[5]!r}")
        intervals.append(interval)
    return intervals


def chrom_names(src: str) -> set[str]:
    """BED names that refer to a MAF source.

    Either the full name, or the part after the species: ``hg38.chr1`` is matched
    by ``hg38.chr1`` and ``chr1``. Versioned accessions are matched in full only.
    """
    names = {src}
    if "." in src and not _ACCESSION.fullmatch(src):
        names.add(src.split(".", 1)[1])
    return names


class IntervalIndex:
    """Find the intervals overlapping a region without scanning all of them."""

    def __init__(self, intervals: Iterable[Interval]):
        grouped: dict[str, list[tuple[int, int, int]]] = {}
        for order, (chrom, start, end, *_) in enumerate(intervals):
            grouped.setdefault(chrom, []).append((start, end, order))
        self._groups = {}
        for chrom, entries in grouped.items():
            entries.sort()
            starts = [e[0] for e in entries]
            max_ends = list(accumulate((e[1] for e in entries), max))
            self._groups[chrom] = (entries, starts, max_ends)

    def overlapping(self, src: str, start: int, end: int) -> list[tuple[int, int, int]]:
        """Overlapping ``(order, start, end)`` entries, in input order."""
        hits = []
        for chrom in chrom_names(src):
            if chrom not in self._groups:
                continue
            entries, starts, max_ends = self._groups[chrom]
            i = bisect_left(starts, end) - 1
            while i >= 0 and max_ends[i] > start:
                s, e, order = entries[i]
                if max(s, start) < min(e, end):
                    hits.append((order, s, e))
                i -= 1
        return sorted(hits)
