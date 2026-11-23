"""Per-block summary measures."""

from __future__ import annotations

from .blocks import GAP, Block


def pairwise_identity(seq1: str, seq2: str) -> float:
    """Case-insensitive identity over columns where not both rows have a gap.

    A base facing a gap counts as a mismatch.
    """
    if len(seq1) != len(seq2):
        raise ValueError("sequences must be of equal length")
    matches = compared = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a == GAP and b == GAP:
            continue
        compared += 1
        matches += a == b
    return matches / compared if compared else 0.0


def reference_identities(block: Block) -> list[float]:
    """Identity of every other row to the reference row."""
    ref = block.reference.text
    return [pairwise_identity(ref, row.text) for row in block.rows[1:]]


def mean_gap_fraction(block: Block) -> float:
    if not block.columns:
        return 0.0
    gaps = sum(row.text.count(GAP) for row in block.rows)
    return gaps / (block.columns * len(block.rows))
