"""Random access to sequences in uncompressed FASTA files."""

from __future__ import annotations

import os
from typing import NamedTuple

from .files import is_gzip
from .mafio import MafError


class Record(NamedTuple):
    length: int
    offset: int  # of the first base in the file
    bases: int  # per line
    width: int  # of a line, including the line end


class Fasta:
    """A FASTA file. Uses its ``.fai`` index, or scans the file once."""

    def __init__(self, path: str):
        self.path = path
        if is_gzip(path):
            raise MafError(f"{path}: compressed FASTA is not supported, please decompress it")
        self._handle = open(path, "rb")
        index = path + ".fai"
        if os.path.exists(index) and os.stat(index).st_mtime_ns >= os.stat(path).st_mtime_ns:
            self.records = self._read_index(index)
        else:
            self.records = self._scan()

    def close(self) -> None:
        self._handle.close()

    def _read_index(self, index: str) -> dict[str, Record]:
        records = {}
        with open(index, encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                try:
                    name, *numbers = line.rstrip("\n").split("\t")
                    records[name] = Record(*map(int, numbers[:4]))
                except (TypeError, ValueError):
                    raise MafError(f"{index}:{number}: malformed index line") from None
        return records

    def _scan(self) -> dict[str, Record]:
        records: dict[str, Record] = {}
        name = None
        length = offset = bases = width = position = 0
        short = False  # a shorter line was seen; only the last one may be short

        def store():
            if name is not None:
                records[name] = Record(length, offset, bases, width)

        for line in self._handle:
            position += len(line)
            if line.startswith(b">"):
                store()
                name = (line[1:].split() or [b""])[0].decode("utf-8")
                if name in records:
                    raise MafError(f"{self.path}: sequence {name!r} occurs twice")
                length, offset, bases, width, short = 0, position, 0, 0, False
                continue
            count = len(line.rstrip())
            if not count:
                short = True  # a blank line, fine at the end of a record
                continue
            if name is None:
                raise MafError(f"{self.path}: not a FASTA file")
            if not bases:
                bases, width = count, len(line)
            if short or count > bases:
                raise MafError(f"{self.path}: lines of {name!r} differ in length")
            short = count < bases
            length += count
        store()
        return records

    def fetch(self, name: str, start: int, end: int) -> str:
        """Upper-case bases ``start:end`` of a sequence, on its plus strand."""
        record = self.records[name]
        start, end = max(0, start), min(end, record.length)
        if start >= end:
            return ""
        first = record.offset + start // record.bases * record.width + start % record.bases
        last = record.offset + end // record.bases * record.width + end % record.bases
        self._handle.seek(first)
        return b"".join(self._handle.read(last - first).split()).decode("utf-8").upper()


FASTA_SUFFIXES = (".fa", ".fasta", ".fna")


class Genomes:
    """Find the sequence that belongs to a MAF source.

    Files given for a species hold its sequences, named either as in the MAF file
    (``hg38.chr1``) or without the species (``chr1``). Other files are searched for
    the full name only.
    """

    def __init__(self) -> None:
        self._by_species: dict[str, list[Fasta]] = {}
        self._shared: list[Fasta] = []

    def add(self, path: str, species: str | None = None) -> None:
        fasta = Fasta(path)
        if species is None:
            self._shared.append(fasta)
        else:
            self._by_species.setdefault(species, []).append(fasta)

    def add_directory(self, directory: str) -> int:
        """Add ``<species>.fa`` (or .fasta, .fna) files; returns their number."""
        count = 0
        for name in sorted(os.listdir(directory)):
            species, suffix = os.path.splitext(name)
            if suffix in FASTA_SUFFIXES:
                self.add(os.path.join(directory, name), species)
                count += 1
        return count

    def files(self) -> list[Fasta]:
        return [f for group in (*self._by_species.values(), self._shared) for f in group]

    def close(self) -> None:
        for fasta in self.files():
            fasta.close()

    def find(self, src: str) -> tuple[Fasta, str] | None:
        """The file and the name of a source, or None."""
        species, _, rest = src.partition(".")
        for fasta in self._by_species.get(species, ()):
            for name in (src, rest):
                if name in fasta.records:
                    return fasta, name
        for fasta in self._shared:
            if src in fasta.records:
                return fasta, src
        return None
