"""In-memory representation of MAF alignment blocks.

Coordinates follow the MAF convention: ``start`` is 0-based, and for rows on the
minus strand it counts from the start of the *reverse-complemented* source.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

GAP = "-"


def count_bases(text: str) -> int:
    """Number of non-gap characters in an aligned sequence."""
    return len(text) - text.count(GAP)


def column_of_base(text: str, index: int) -> int:
    """Column holding the ``index``-th (0-based) non-gap character."""
    seen = 0
    for column, char in enumerate(text):
        if char != GAP:
            if seen == index:
                return column
            seen += 1
    raise IndexError(f"base {index} out of range ({seen} bases in row)")


@dataclass(slots=True)
class Row:
    """One aligned sequence: an ``s`` line plus its optional ``q``/``i`` lines."""

    src: str
    start: int
    size: int
    strand: str
    src_size: int
    text: str
    quality: str | None = None
    info: str | None = None

    @property
    def end(self) -> int:
        """Exclusive end, on the row's own strand."""
        return self.start + self.size

    @property
    def species(self) -> str:
        return self.src.split(".", 1)[0]

    def plus_interval(self) -> tuple[int, int]:
        """Half-open interval on the plus strand, as used by BED."""
        if self.strand == "+":
            return self.start, self.end
        return self.src_size - self.end, self.src_size - self.start

    def to_strand(self, start: int, end: int) -> tuple[int, int]:
        """Convert a plus-strand interval to this row's strand."""
        if self.strand == "+":
            return start, end
        return self.src_size - end, self.src_size - start

    def sliced(self, first: int, last: int) -> Row:
        """Columns ``first:last``, with start and size adjusted to match."""
        text = self.text[first:last]
        return replace(
            self,
            start=self.start + count_bases(self.text[:first]),
            size=count_bases(text),
            text=text,
            quality=None if self.quality is None else self.quality[first:last],
            info=None,
        )


@dataclass(slots=True)
class Block:
    """One alignment block. The first row is the reference."""

    rows: list[Row] = field(default_factory=list)
    attrs: dict[str, str | None] = field(default_factory=dict)
    empties: list[str] = field(default_factory=list)  # fields of raw `e` lines

    @property
    def reference(self) -> Row:
        return self.rows[0]

    @property
    def columns(self) -> int:
        return len(self.rows[0].text) if self.rows else 0

    def derived(self, rows: list[Row]) -> Block:
        """A block with edited rows.

        The score and the ``i``/``e`` context lines describe the original block
        and no longer hold, so they are not carried over.
        """
        attrs = {k: v for k, v in self.attrs.items() if k != "score"}
        return Block([replace(row, info=None) for row in rows], attrs)

    def sliced(self, first: int, last: int) -> Block | None:
        """Columns ``first:last``; rows without any base there are dropped.

        Returns None if the reference has no base in the slice, because such a
        block has no reference position.
        """
        if first <= 0 and last >= self.columns:
            return self
        rows = [row.sliced(first, last) for row in self.rows]
        if not rows[0].size:
            return None
        return self.derived([row for row in rows if row.size])

    def without_gap_columns(self) -> Block:
        """Remove columns in which every row has a gap."""
        texts = [row.text for row in self.rows]
        keep = [any(c != GAP for c in column) for column in zip(*texts)]
        if all(keep):
            return self

        def squeeze(text: str) -> str:
            return "".join(c for c, kept in zip(text, keep) if kept)

        return self.derived(
            [
                replace(
                    row,
                    text=squeeze(row.text),
                    quality=None if row.quality is None else squeeze(row.quality),
                )
                for row in self.rows
            ]
        )
