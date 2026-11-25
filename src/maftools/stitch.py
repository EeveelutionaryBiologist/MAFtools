"""One alignment per feature, put together from all blocks that overlap it."""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from .bed import Interval, IntervalIndex
from .blocks import GAP, Block, reverse_complement
from .operations import reference_columns


@dataclass
class Stitched:
    """Aligned texts by species, the reference first, in the feature's direction."""

    name: str
    src: str
    start: int
    end: int
    strand: str
    texts: dict[str, str]
    covered: int  # reference bases in alignment blocks

    @property
    def location(self) -> str:
        return f"{self.src}:{self.start}-{self.end}({self.strand})"


def feature_name(interval: Interval) -> str:
    """The BED name, or the position if there is none."""
    if interval.name not in ("", "."):
        return re.sub(r"\s+", "_", interval.name)
    return f"{interval.chrom}_{interval.start}_{interval.end}_{interval.strand}"


def _join(parts: list[Block], start: int, end: int, stats: Counter) -> tuple[dict[str, str], int]:
    """Join plus-strand parts. Reference bases without alignment become N."""
    reference = parts[0].reference.species
    texts: dict[str, list[str]] = {}
    lengths: Counter = Counter()
    total = covered = 0
    position = start

    def append(species: str, text: str) -> None:
        texts.setdefault(species, []).extend((GAP * (total - lengths[species]), text))
        lengths[species] = total + len(text)

    def fill(until: int) -> None:
        nonlocal total
        append(reference, "N" * (until - position))
        total += until - position

    for part in sorted(parts, key=lambda p: p.reference.start):  # stable: file order
        if part.reference.end <= position:
            stats["overlaps"] += 1
            continue
        if part.reference.start < position:
            stats["overlaps"] += 1
            part = part.sliced(*reference_columns(part.reference, position, part.reference.end))
        fill(part.reference.start)
        used = set()
        for row in part.rows:
            if row.species in used:
                stats["duplicates"] += 1
                continue
            used.add(row.species)
            append(row.species, row.text)
        total += part.columns
        covered += part.reference.size
        position = part.reference.end
    fill(end)
    return (
        {s: "".join(t) + GAP * (total - lengths[s]) for s, t in texts.items()},
        covered,
    )


def stitch(
    blocks: Iterable[Block], intervals: Iterable[Interval], stats: Counter | None = None
) -> Iterator[Stitched]:
    """One alignment for every interval with an alignment, in input order.

    Without gaps, the reference has exactly the length of the interval. Where
    blocks overlap on the reference, the one that starts first is used. If a
    species occurs more than once in a block, its first row is used. Intervals on
    the minus strand are reverse-complemented.
    """
    stats = Counter() if stats is None else stats
    intervals = list(intervals)
    stats["intervals"] = len(intervals)
    index = IntervalIndex(intervals)
    found: dict[tuple[int, str], list[Block]] = {}
    for block in blocks:
        ref = block.reference
        for order, start, end in index.overlapping(ref.src, *ref.plus_interval()):
            part = block.sliced(*reference_columns(ref, *ref.to_strand(start, end)))
            if part is None:
                continue
            if part.reference.strand == "-":
                part = part.reverse_complement()
            found.setdefault((order, ref.src), []).append(part)

    for order, src in sorted(found):
        interval = intervals[order]
        # a feature may reach beyond the end of its sequence
        end = min(interval.end, found[order, src][0].reference.src_size)
        texts, covered = _join(found[order, src], interval.start, end, stats)
        if interval.strand == "-":
            texts = {s: reverse_complement(t) for s, t in texts.items()}
        stats["partial"] += covered < end - interval.start
        yield Stitched(
            feature_name(interval), src, interval.start, end, interval.strand, texts, covered
        )
    stats["matched"] = len({order for order, _ in found})


def file_names(names: Iterable[str]) -> Iterator[str]:
    """Distinct names that are safe as file names."""
    taken: set[str] = set()
    for name in names:
        # no leading dot: no hidden files, and no way out of the directory
        base = re.sub(r"[^A-Za-z0-9._+-]", "_", name).lstrip(".")[:200] or "_"
        candidate, n = base, 1
        # not all file systems distinguish case
        while candidate.lower() in taken:
            n += 1
            candidate = f"{base}_{n}"
        taken.add(candidate.lower())
        yield candidate
