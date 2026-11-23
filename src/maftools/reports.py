"""Summaries and conversions that do not produce MAF."""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from .blocks import Block
from .metrics import mean_gap_fraction, reference_identities

BLOCK_COLUMNS = (
    "refseq", "start", "size", "length", "seqs", "mean_pairwise_identity", "gaps_fraction",
)


def describe(blocks: Iterable[Block]) -> dict[str, int | float]:
    """File-level statistics. Averages are 0 for a file without blocks."""
    species, sources = set(), set()
    n_blocks = columns = bases = rows = 0
    for block in blocks:
        n_blocks += 1
        columns += block.columns
        bases += block.reference.size
        rows += len(block.rows)
        sources.update(row.src for row in block.rows)
        species.update(row.species for row in block.rows)
    return {
        "Identifiers": len(species),
        "Sequences": len(sources),
        "Number of blocks": n_blocks,
        "Columns": columns,
        "Average sequences per block": round(rows / n_blocks, 2) if n_blocks else 0,
        "Aligned reference nucleotides": bases,
        "Average block length": round(columns / n_blocks, 2) if n_blocks else 0,
    }


def block_table(blocks: Iterable[Block]) -> Iterator[str]:
    """Tab-separated lines, header first, one per block."""
    yield "\t".join(BLOCK_COLUMNS)
    for block in blocks:
        ref = block.reference
        identities = reference_identities(block)
        mean = sum(identities) / len(identities) if identities else 1.0
        yield (
            f"{ref.src}\t{ref.start}\t{ref.size}\t{block.columns}\t{len(block.rows)}\t"
            f"{mean:.4f}\t{mean_gap_fraction(block):.4f}"
        )


def bed_lines(blocks: Iterable[Block]) -> Iterator[str]:
    """Plus-strand reference interval of every block."""
    for block in blocks:
        start, end = block.reference.plus_interval()
        yield f"{block.reference.src}\t{start}\t{end}"
