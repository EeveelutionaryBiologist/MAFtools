"""Streaming transformations: each takes blocks and yields blocks.

``stats`` is an optional Counter that collects what happened, for reporting.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, replace

from .bed import Interval, IntervalIndex
from .blocks import GAP, Block, Row, column_of_base
from .metrics import mean_gap_fraction, pairwise_identity, reference_identities


# --- select -----------------------------------------------------------------

def reference_columns(row: Row, start: int, end: int) -> tuple[int, int]:
    """Columns covering bases ``start:end`` (the row's strand) of a row.

    Insertions relative to the row belong to the base before them, and those at
    the very start of a block to its first base, so that adjacent intervals
    partition the columns of a block.
    """
    first = 0 if start <= row.start else column_of_base(row.text, start - row.start)
    if end >= row.end:
        return first, len(row.text)
    return first, column_of_base(row.text, end - row.start)


def widen(intervals: Iterable[Interval], left: int, right: int) -> list[Interval]:
    """Add flanks on the plus strand."""
    return [Interval(c, max(0, s - left), e + right) for c, s, e in intervals]


def select_blocks(
    blocks: Iterable[Block],
    intervals: Iterable[Interval],
    left: int = 0,
    right: int = 0,
    stats: Counter | None = None,
) -> Iterator[Block]:
    """Extract the parts of blocks that overlap intervals on the reference.

    Flanks extend intervals on the plus strand; the BED strand is not used.
    """
    stats = Counter() if stats is None else stats
    widened = widen(intervals, left, right)
    index = IntervalIndex(widened)
    matched: set[int] = set()
    stats["intervals"] = len(widened)
    for block in blocks:
        ref = block.reference
        for order, start, end in index.overlapping(ref.src, *ref.plus_interval()):
            first, last = reference_columns(ref, *ref.to_strand(start, end))
            part = block.sliced(first, last)
            if part is not None:
                matched.add(order)
                stats["matched"] = len(matched)
                yield part


# --- filter -----------------------------------------------------------------

@dataclass(frozen=True)
class FilterCriteria:
    max_gaps: float = 0.25
    min_similarity: float = 0.5
    max_similarity: float = 1.0
    min_length: int = 40
    min_seqs: int = 3

    def accepts(self, block: Block) -> bool:
        """Every row must lie within the similarity bounds, not just the mean."""
        return (
            len(block.rows) >= self.min_seqs
            and block.columns >= self.min_length
            and mean_gap_fraction(block) <= self.max_gaps
            and all(
                self.min_similarity <= identity <= self.max_similarity
                for identity in reference_identities(block)
            )
        )


def filter_blocks(
    blocks: Iterable[Block], criteria: FilterCriteria, stats: Counter | None = None
) -> Iterator[Block]:
    stats = Counter() if stats is None else stats
    for block in blocks:
        stats["read"] += 1
        if criteria.accepts(block):
            stats["kept"] += 1
            yield block


# --- species ----------------------------------------------------------------

@dataclass(frozen=True)
class SpeciesOptions:
    """Names are species (``hg38``) or full sequence names (``hg38.chr1``)."""

    names: tuple[str, ...]
    exclude: bool = False
    require_all: bool = False
    reorder: bool = False  # rows in the order of names; the first is the reference
    min_seqs: int = 1

    def rank(self, row: Row) -> int | None:
        """Position of the first name that matches a row."""
        for i, name in enumerate(self.names):
            if name in (row.species, row.src):
                return i
        return None


def species_blocks(
    blocks: Iterable[Block],
    options: SpeciesOptions,
    stats: Counter | None = None,
    seen: set[str] | None = None,
) -> Iterator[Block]:
    """Keep, drop or reorder rows. Columns left without any base are removed.

    With ``reorder``, the first name becomes the reference: blocks without it are
    dropped and blocks with it on the minus strand are reverse-complemented.
    ``seen`` collects the names that matched any row.
    """
    stats = Counter() if stats is None else stats
    seen = set() if seen is None else seen
    for block in blocks:
        stats["read"] += 1
        ranks = [options.rank(row) for row in block.rows]
        found = {rank for rank in ranks if rank is not None}
        seen.update(options.names[rank] for rank in found)
        if options.require_all and len(found) < len(options.names):
            continue
        if options.reorder and 0 not in found:
            continue
        kept = [i for i, rank in enumerate(ranks) if (rank is None) == options.exclude]
        if options.reorder:
            kept.sort(key=lambda i: ranks[i])
        if not kept or len(kept) < options.min_seqs:
            continue
        if kept != list(range(len(ranks))):
            stats["changed_reference"] += kept[0] != 0
            block = block.derived([block.rows[i] for i in kept]).without_gap_columns()
        if options.reorder and block.reference.strand == "-":
            block = block.reverse_complement()
        if block.reference.size:
            stats["kept"] += 1
            yield block


# --- dedup ------------------------------------------------------------------

KEEP_BEST, KEEP_FIRST, KEEP_NONE = "best", "first", "none"


def _copies_to_keep(block: Block, copies: list[int], keep: str) -> list[int]:
    """Which of several rows of one species stay. The reference always wins."""
    if keep == KEEP_NONE:
        return []
    if keep == KEEP_FIRST or copies[0] == 0:
        return copies[:1]
    ref = block.reference.text
    # max() returns the first of equally good copies
    return [max(copies, key=lambda i: pairwise_identity(ref, block.rows[i].text))]


def dedup_blocks(
    blocks: Iterable[Block],
    keep: str = KEEP_BEST,
    min_seqs: int = 1,
    stats: Counter | None = None,
) -> Iterator[Block]:
    """Reduce blocks to at most one row per species.

    ``best`` keeps the copy most similar to the reference, ``first`` the first one,
    and ``none`` drops all copies. With ``none``, blocks in which the species of the
    reference is duplicated are dropped. Columns left without any base are removed.
    """
    if keep not in (KEEP_BEST, KEEP_FIRST, KEEP_NONE):
        raise ValueError(f"unknown value for keep: {keep!r}")
    stats = Counter() if stats is None else stats
    for block in blocks:
        stats["read"] += 1
        species: dict[str, list[int]] = {}
        for i, row in enumerate(block.rows):
            species.setdefault(row.species, []).append(i)
        kept: list[int] = []
        for copies in species.values():
            kept += copies if len(copies) == 1 else _copies_to_keep(block, copies, keep)
        kept.sort()
        removed = len(block.rows) - len(kept)
        stats["affected"] += removed > 0
        if 0 not in kept or len(kept) < min_seqs:
            stats["dropped"] += 1
            continue
        if removed:
            stats["removed"] += removed
            block = block.derived([block.rows[i] for i in kept]).without_gap_columns()
        stats["kept"] += 1
        yield block


# --- window -----------------------------------------------------------------

def window_starts(columns: int, length: int, step: int) -> list[int]:
    """Window starts covering all columns; the last window is shifted to fit."""
    if columns <= length:
        return [0]
    starts = list(range(0, columns - length + 1, step))
    if starts[-1] + length < columns:
        starts.append(columns - length)
    return starts


def window_blocks(
    blocks: Iterable[Block], length: int, step: int, stats: Counter | None = None
) -> Iterator[Block]:
    """Cut blocks into overlapping windows of ``length`` columns.

    Windows in which the reference has no base are skipped.
    """
    if not 0 < step <= length:
        raise ValueError("step must be positive and not exceed the window length")
    stats = Counter() if stats is None else stats
    for block in blocks:
        stats["read"] += 1
        for start in window_starts(block.columns, length, step):
            window = block.sliced(start, start + length)
            if window is None:
                stats["skipped"] += 1
            else:
                stats["windows"] += 1
                yield window


# --- merge ------------------------------------------------------------------

@dataclass(frozen=True)
class MergeOptions:
    species_consensus: float = 0.75
    max_distance: int = 0
    max_length: int = 1000
    max_gaps: float = 0.5
    min_seqs: int = 2
    keep_reference: bool = True  # exempt the first row from max_gaps


def _distance(left: Row, right: Row) -> int:
    return right.start - left.end


def _joinable(left: Row, right: Row, max_distance: int) -> bool:
    return (
        left.src == right.src
        and left.strand == right.strand
        and 0 <= _distance(left, right) <= max_distance
    )


def _pair_rows(a: Block, b: Block, max_distance: int) -> list[tuple[Row, Row]]:
    """Match rows by name and position, never by their order in the block.

    The references always form the first pair. If a source occurs more than once,
    the closest copies are paired first. No row is used twice.
    """
    candidates = sorted(
        (_distance(left, right), i, j)
        for i, left in enumerate(a.rows[1:], 1)
        for j, right in enumerate(b.rows[1:], 1)
        if _joinable(left, right, max_distance)
    )
    partner = {0: 0}
    used = {0}
    for _, i, j in candidates:
        if i not in partner and j not in used:
            partner[i] = j
            used.add(j)
    return [(a.rows[i], b.rows[j]) for i, j in sorted(partner.items())]


def merge_pair(a: Block, b: Block, options: MergeOptions) -> Block | None:
    """Join two consecutive blocks, or return None if they must stay apart.

    Rows that are not continued in the second block are dropped. Unaligned bases
    between the blocks become ``N``, padded with gaps to equal length.
    """
    if not _joinable(a.reference, b.reference, options.max_distance):
        return None
    ids_a = {row.src for row in a.rows}
    ids_b = {row.src for row in b.rows}
    consensus = 2 * len(ids_a & ids_b) / (len(ids_a) + len(ids_b))
    if consensus < options.species_consensus:
        return None
    pairs = _pair_rows(a, b, options.max_distance)
    if len(pairs) < options.min_seqs:
        return None
    bridge = max(_distance(left, right) for left, right in pairs)
    if a.columns + bridge + b.columns > options.max_length:
        return None
    rows = []
    for left, right in pairs:
        gap = _distance(left, right)
        rows.append(
            replace(
                left,
                size=left.size + gap + right.size,
                text=left.text + "N" * gap + GAP * (bridge - gap) + right.text,
                quality=None,
            )
        )
    return a.derived(rows)


def _finish(block: Block, options: MergeOptions) -> Block | None:
    """Drop gap-only columns and gappy rows before a block is written."""
    block = block.without_gap_columns()
    if not block.columns:
        return None
    rows = [
        row
        for i, row in enumerate(block.rows)
        if (i == 0 and options.keep_reference)
        or (row.size and row.text.count(GAP) / block.columns <= options.max_gaps)
    ]
    if len(rows) < options.min_seqs or not rows:
        return None
    if len(rows) < len(block.rows):
        block = block.derived(rows).without_gap_columns()
    return block


def merge_blocks(
    blocks: Iterable[Block], options: MergeOptions, stats: Counter | None = None
) -> Iterator[Block]:
    """Merge runs of consecutive blocks."""
    stats = Counter() if stats is None else stats

    def finish(block: Block) -> Iterator[Block]:
        done = _finish(block, options)
        if done is None:
            stats["dropped"] += 1
        else:
            stats["written"] += 1
            yield done

    current: Block | None = None
    for block in blocks:
        stats["read"] += 1
        if current is None:
            current = block
            continue
        merged = merge_pair(current, block, options)
        if merged is None:
            yield from finish(current)
            current = block
        else:
            stats["merges"] += 1
            current = merged
    if current is not None:
        yield from finish(current)
