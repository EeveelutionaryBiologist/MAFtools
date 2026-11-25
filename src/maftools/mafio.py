"""Streaming MAF reader and writer.

Biopython's MAF support drops block scores and ``i``/``e``/``q`` lines, so the
format is handled here directly. One block is held in memory at a time.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from itertools import chain
from typing import TextIO

from .blocks import Block, Row, count_bases

DEFAULT_HEADER = "##maf version=1"


class MafError(ValueError):
    """Malformed input."""

    def __init__(self, message: str, line: int = 0):
        super().__init__(message)
        self.line = line
        self.reason = message


class MafReader:
    """Iterate over the blocks of a MAF file.

    ``header`` holds the comment lines that precede the first block. Lines of
    unknown type are ignored, as the MAF specification asks.

    If ``on_error`` is set, it receives the errors of malformed lines, which are
    then skipped. Otherwise the first error is raised.
    """

    def __init__(self, handle: Iterable[str], name: str = "<maf>"):
        self.name = name
        self.on_error: Callable[[MafError], None] | None = None
        self.block_line = 0  # line number of the block yielded last
        self.header: list[str] = []
        self._lines = enumerate(handle, 1)
        self._first: list[tuple[int, str]] = []
        for number, line in self._lines:
            if line.startswith("#"):
                self.header.append(line.rstrip())
            elif line.strip() and not line.startswith("track"):
                self._first = [(number, line)]
                break

    def _error(self, number: int, message: str) -> MafError:
        error = MafError(f"{self.name}:{number}: {message}", number)
        error.reason = message
        return error

    def __iter__(self) -> Iterator[Block]:
        block: Block | None = None
        line_of_block = 0
        for number, line in chain(self._first, self._lines):
            fields = line.split()
            if not fields or fields[0] == "a":
                if block is not None and block.rows:
                    self.block_line = line_of_block
                    yield block
                block = None
                if fields:
                    line_of_block = number
                    attrs = (f.partition("=") for f in fields[1:])
                    block = Block(attrs={k: v if eq else None for k, eq, v in attrs})
            elif fields[0].startswith("#"):
                continue
            elif block is None:
                self._handle(self._error(number, "expected an alignment block ('a' line)"))
            elif fields[0] == "s":
                try:
                    block.rows.append(self._row(number, fields, block))
                except MafError as error:
                    self._handle(error)
            elif fields[0] in ("i", "q"):
                try:
                    self._annotate(number, fields, block)
                except MafError as error:
                    self._handle(error)
            elif fields[0] == "e":
                block.empties.append(" ".join(fields[1:]))
        if block is not None and block.rows:
            self.block_line = line_of_block
            yield block

    def _handle(self, error: MafError) -> None:
        if self.on_error is None:
            raise error
        self.on_error(error)

    def _row(self, number: int, fields: list[str], block: Block) -> Row:
        if len(fields) != 7:
            raise self._error(number, f"'s' line has {len(fields)} fields, expected 7")
        _, src, start, size, strand, src_size, text = fields
        try:
            row = Row(src, int(start), int(size), strand, int(src_size), text)
        except ValueError:
            raise self._error(number, "start, size and srcSize must be integers") from None
        if strand not in ("+", "-"):
            raise self._error(number, f"invalid strand {strand!r}")
        if row.size != count_bases(text):
            raise self._error(
                number, f"size is {size} but the sequence has {count_bases(text)} bases"
            )
        if block.rows and len(text) != block.columns:
            raise self._error(number, "rows of a block differ in length")
        return row

    def _annotate(self, number: int, fields: list[str], block: Block) -> None:
        if len(fields) < 3 or not block.rows or block.rows[-1].src != fields[1]:
            raise self._error(number, f"'{fields[0]}' line does not follow its 's' line")
        row = block.rows[-1]
        if fields[0] == "i":
            row.info = " ".join(fields[2:])
        elif len(fields[2]) != len(row.text):
            raise self._error(number, "'q' line and sequence differ in length")
        else:
            row.quality = fields[2]


def format_block(block: Block) -> str:
    """Render a block as MAF text, including the trailing blank line."""
    attrs = (k if v is None else f"{k}={v}" for k, v in block.attrs.items())
    lines = [" ".join(["a", *attrs])]
    rows = block.rows
    src_w = max(len(r.src) for r in rows)
    start_w = max(len(str(r.start)) for r in rows)
    size_w = max(len(str(r.size)) for r in rows)
    total_w = max(len(str(r.src_size)) for r in rows)
    for r in rows:
        lines.append(
            f"s {r.src:<{src_w}} {r.start:>{start_w}} {r.size:>{size_w}} "
            f"{r.strand} {r.src_size:>{total_w}} {r.text}"
        )
        if r.quality is not None:
            lines.append(f"q {r.src} {r.quality}")
        if r.info is not None:
            lines.append(f"i {r.src} {r.info}")
    lines.extend(f"e {e}" for e in block.empties)
    return "\n".join(lines) + "\n\n"


def write_maf(handle: TextIO, blocks: Iterable[Block], header: Iterable[str] = ()) -> int:
    """Write a complete MAF file and return the number of blocks."""
    header = list(header)
    if not header or not header[0].startswith("##maf"):
        header.insert(0, DEFAULT_HEADER)
    handle.write("\n".join(header) + "\n\n")
    written = 0
    for block in blocks:
        handle.write(format_block(block))
        written += 1
    return written
