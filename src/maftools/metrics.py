"""Per-block summary measures."""

from __future__ import annotations

from .blocks import GAP, Block


class Encoded:
    """An aligned sequence as integers with one byte per column.

    Comparing integers is much faster than comparing strings column by column.
    """

    __slots__ = ("columns", "bases", "gaps")
    _IS_GAP = bytes(char == ord(GAP) for char in range(256))

    def __init__(self, text: str):
        data = text.upper().encode("latin-1", "replace")
        self.columns = len(data)
        self.bases = int.from_bytes(data, "little")
        self.gaps = int.from_bytes(data.translate(self._IS_GAP), "little")

    def identity(self, other: Encoded) -> float:
        n = self.columns
        if n != other.columns:
            raise ValueError("sequences must be of equal length")
        shared_gaps = (self.gaps & other.gaps).to_bytes(n, "little").count(1)
        matches = (self.bases ^ other.bases).to_bytes(n, "little").count(0) - shared_gaps
        return matches / (n - shared_gaps) if n > shared_gaps else 0.0


def pairwise_identity(seq1: str, seq2: str) -> float:
    """Case-insensitive identity over columns where not both rows have a gap.

    A base facing a gap counts as a mismatch.
    """
    return Encoded(seq1).identity(Encoded(seq2))


def reference_identities(block: Block) -> list[float]:
    """Identity of every other row to the reference row."""
    ref = Encoded(block.reference.text)
    return [ref.identity(Encoded(row.text)) for row in block.rows[1:]]


def mean_gap_fraction(block: Block) -> float:
    if not block.columns:
        return 0.0
    gaps = sum(row.text.count(GAP) for row in block.rows)
    return gaps / (block.columns * len(block.rows))
