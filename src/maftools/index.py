"""On-disk index of reference intervals, for random access to large MAF files.

The index is an SQLite database next to the MAF file. It records the size and
modification time of the MAF file, so that an outdated index is recognised.
"""

from __future__ import annotations

import os
import sqlite3
from collections.abc import Iterable, Iterator
from contextlib import closing
from pathlib import Path

from .bed import Interval, chrom_names
from .blocks import Block
from .files import STDIO, is_gzip
from .mafio import MafError, MafReader

SUFFIX = ".mfi"
VERSION = 1

_SCHEMA = """
CREATE TABLE meta (version INTEGER, size INTEGER, mtime INTEGER);
CREATE TABLE sources (id INTEGER PRIMARY KEY, src TEXT UNIQUE, longest INTEGER);
CREATE TABLE blocks (source INTEGER, start INTEGER, end INTEGER, offset INTEGER);
"""


def index_path(maf: str) -> str:
    return maf + SUFFIX


def _fingerprint(maf: str) -> tuple[int, int]:
    info = os.stat(maf)
    return info.st_size, info.st_mtime_ns


class _Lines:
    """Decoded lines of a binary file that remember where 'a' lines start."""

    def __init__(self, handle):
        self._handle = handle
        self.offsets: dict[int, int] = {}

    def __iter__(self) -> Iterator[str]:
        offset = 0
        for number, line in enumerate(self._handle, 1):
            if line.split(None, 1)[:1] == [b"a"]:
                self.offsets[number] = offset
            offset += len(line)
            yield line.decode("utf-8")


def build_index(maf: str) -> int:
    """Index the reference interval of every block; returns the block count."""
    if maf == STDIO or is_gzip(maf):
        raise MafError("only uncompressed files can be indexed")
    target = index_path(maf)
    temporary = f"{target}.{os.getpid()}.tmp"
    fingerprint = _fingerprint(maf)
    sources: dict[str, int] = {}
    longest: dict[int, int] = {}

    def entries(reader: MafReader, lines: _Lines):
        for block in reader:
            offset = lines.offsets.pop(reader.block_line)
            ref = block.reference
            start, end = ref.plus_interval()
            source = sources.setdefault(ref.src, len(sources))
            longest[source] = max(longest.get(source, 0), end - start)
            yield source, start, end, offset

    try:
        with open(maf, "rb") as handle, closing(sqlite3.connect(temporary)) as db:
            lines = _Lines(handle)
            db.executescript(_SCHEMA)
            cursor = db.executemany(
                "INSERT INTO blocks VALUES (?, ?, ?, ?)", entries(MafReader(lines, maf), lines)
            )
            count = max(cursor.rowcount, 0)
            db.executemany(
                "INSERT INTO sources VALUES (?, ?, ?)",
                [(i, src, longest[i]) for src, i in sources.items()],
            )
            db.execute("CREATE INDEX by_start ON blocks (source, start, end, offset)")
            db.execute("INSERT INTO meta VALUES (?, ?, ?)", (VERSION, *fingerprint))
            db.commit()
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)
    return count


class MafIndex:
    """Read blocks overlapping intervals, without reading the rest of the file."""

    def __init__(self, maf: str):
        self.maf = maf
        uri = Path(index_path(maf)).resolve().as_uri() + "?mode=ro"
        self._db = sqlite3.connect(uri, uri=True)

    def close(self) -> None:
        self._db.close()

    def is_current(self) -> bool:
        try:
            row = self._db.execute("SELECT version, size, mtime FROM meta").fetchone()
        except sqlite3.DatabaseError:
            return False
        return row == (VERSION, *_fingerprint(self.maf))

    def offsets(self, intervals: Iterable[Interval]) -> list[int]:
        """File offsets of overlapping blocks, in file order."""
        by_name: dict[str, list[tuple[int, int]]] = {}
        for source, src, longest in self._db.execute("SELECT id, src, longest FROM sources"):
            for name in chrom_names(src):
                by_name.setdefault(name, []).append((source, longest))
        found: set[int] = set()
        for chrom, start, end in intervals:
            for source, longest in by_name.get(chrom, ()):
                found.update(
                    offset
                    for (offset,) in self._db.execute(
                        "SELECT offset FROM blocks WHERE source = ? AND start < ? "
                        "AND start > ? AND end > ?",
                        (source, end, start - longest, start),
                    )
                )
        return sorted(found)

    def blocks(self, intervals: Iterable[Interval]) -> Iterator[Block]:
        with open(self.maf, "rb") as handle:
            for offset in self.offsets(intervals):
                handle.seek(offset)
                lines = (line.decode("utf-8") for line in handle)
                yield next(iter(MafReader(lines, f"{self.maf} (byte {offset})")))
