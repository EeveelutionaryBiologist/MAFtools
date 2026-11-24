"""Opening inputs and outputs: plain files, gzip and standard streams."""

from __future__ import annotations

import gzip
import io
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from typing import TextIO

STDIO = "-"
_GZIP_MAGIC = b"\x1f\x8b"


def is_gzip(path: str) -> bool:
    with open(path, "rb") as handle:
        return handle.read(2) == _GZIP_MAGIC


@contextmanager
def open_input(path: str) -> Iterator[TextIO]:
    """Read text from a file, or from standard input for ``-``.

    Gzip (including bgzip) is recognised by content, not by file name.
    """
    raw = sys.stdin.buffer if path == STDIO else open(path, "rb")
    try:
        data = gzip.GzipFile(fileobj=raw) if raw.peek(2)[:2] == _GZIP_MAGIC else raw
        text = io.TextIOWrapper(data, encoding="utf-8")
        yield text
        text.detach()
    finally:
        if path != STDIO:
            raw.close()


@contextmanager
def open_output(path: str | None) -> Iterator[TextIO]:
    """Write text to a file, or to standard output for None or ``-``.

    Files named ``*.gz`` are compressed.
    """
    if path is None or path == STDIO:
        yield sys.stdout
    elif path.endswith(".gz"):
        with open(path, "wb") as raw:
            # no time stamp and file name, so that equal content gives equal files
            with gzip.GzipFile("", "wb", compresslevel=6, fileobj=raw, mtime=0) as packed:
                with io.TextIOWrapper(packed, encoding="utf-8") as text:
                    yield text
    else:
        with open(path, "w", encoding="utf-8") as text:
            yield text
