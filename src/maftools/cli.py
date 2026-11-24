"""Command line interface."""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import zlib
from collections import Counter
from collections.abc import Iterable, Iterator, Sequence
from contextlib import contextmanager
from typing import TextIO

from . import __version__
from .bed import read_bed
from .blocks import Block
from .files import STDIO, open_input, open_output
from .index import MafIndex, build_index, index_path
from .mafio import MafError, MafReader, write_maf
from .operations import (
    FilterCriteria, MergeOptions, SpeciesOptions, filter_blocks, merge_blocks, select_blocks,
    species_blocks, widen, window_blocks,
)
from .reports import bed_lines, block_table, describe, fasta_blocks, fasta_concatenated


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def _get_help_string(self, action):
        if action.default in (None, False, True, []):
            return action.help
        return super()._get_help_string(action)


def _fraction(text: str) -> float:
    value = float(text)
    if not 0 <= value <= 1:
        raise argparse.ArgumentTypeError(f"{text} is not between 0 and 1")
    return value


def _at_least(minimum: int):
    def convert(text: str) -> int:
        value = int(text)
        if value < minimum:
            raise argparse.ArgumentTypeError(f"{text} is smaller than {minimum}")
        return value

    convert.__name__ = "integer"
    return convert


@contextmanager
def _output(args: argparse.Namespace) -> Iterator[TextIO]:
    """The -o file, or stdout. Refuses to overwrite a file that is being read."""
    if args.output not in (None, STDIO) and os.path.exists(args.output):
        for path in (args.maf, getattr(args, "bed", None), getattr(args, "list", None)):
            if path not in (None, STDIO) and os.path.samefile(path, args.output):
                raise MafError(f"{args.output} is also an input file")
    with open_output(args.output) as handle:
        yield handle


def _note(message: str) -> None:
    print(message, file=sys.stderr)


def _write(args: argparse.Namespace, blocks: Iterable[Block], reader: MafReader) -> int:
    with _output(args) as out:
        return write_maf(out, blocks, reader.header)


def _write_lines(args: argparse.Namespace, lines: Iterable[str]) -> None:
    with _output(args) as out:
        out.writelines(line + "\n" for line in lines)


def _describe(args: argparse.Namespace, reader: MafReader) -> None:
    for key, value in describe(reader).items():
        print(f"{key}: {value}")


def _blockdescribe(args: argparse.Namespace, reader: MafReader) -> None:
    _write_lines(args, block_table(reader))


def _to_bed(args: argparse.Namespace, reader: MafReader) -> None:
    _write_lines(args, bed_lines(reader))


def _to_fasta(args: argparse.Namespace, reader: MafReader) -> None:
    if not args.concatenate:
        _write_lines(args, fasta_blocks(reader, args.ungapped))
        return
    stats: Counter = Counter()
    _write_lines(args, fasta_concatenated(reader, stats))
    if stats["duplicates"]:
        _note(f"{stats['duplicates']} additional rows of a species in the same block not used")


def _index(args: argparse.Namespace, reader: MafReader) -> None:
    _note(f"{build_index(args.maf)} blocks indexed in {index_path(args.maf)}")


def _open_index(args: argparse.Namespace) -> MafIndex | None:
    """The index of the input, if there is one that fits."""
    if args.maf == STDIO or not os.path.exists(index_path(args.maf)):
        return None
    index = MafIndex(args.maf)
    if index.is_current():
        return index
    index.close()
    _note(f"warning: {index_path(args.maf)} is out of date and not used")
    return None


def _select(args: argparse.Namespace, reader: MafReader) -> None:
    with open_input(args.bed) as handle:
        intervals = read_bed(handle, args.bed)
    stats: Counter = Counter()
    left, right = args.left_flank, args.right_flank
    blocks: Iterable[Block] = reader
    index = _open_index(args)
    if index is not None:
        blocks = index.blocks(widen(intervals, left, right))
    written = _write(args, select_blocks(blocks, intervals, left, right, stats), reader)
    if index is not None:
        index.close()
    _note(f"{written} blocks written; {stats['matched']} of {stats['intervals']} intervals matched")


def _species(args: argparse.Namespace, reader: MafReader) -> None:
    names = [name for item in args.species for name in item.split(",") if name]
    if args.list is not None:
        with open_input(args.list) as handle:
            names += [f[0] for f in map(str.split, handle) if f and not f[0].startswith("#")]
    if not names:
        raise MafError("no species given, use --species or --list")
    names = list(dict.fromkeys(names))
    options = SpeciesOptions(
        tuple(names), args.exclude, args.require_all, args.reorder, args.min_seqs
    )
    stats: Counter = Counter()
    seen: set[str] = set()
    _write(args, species_blocks(reader, options, stats, seen), reader)
    _note(f"{stats['kept']} of {stats['read']} blocks kept")
    if stats["changed_reference"] and not args.reorder:
        _note(f"warning: {stats['changed_reference']} blocks lost their reference; "
              "the next sequence takes its place")
    missing = [name for name in names if name not in seen]
    if missing:
        _note(f"warning: not found: {', '.join(missing)}")


def _filter(args: argparse.Namespace, reader: MafReader) -> None:
    criteria = FilterCriteria(
        args.max_gaps, args.min_similarity, args.max_similarity, args.min_length, args.min_seqs
    )
    stats: Counter = Counter()
    _write(args, filter_blocks(reader, criteria, stats), reader)
    _note(f"{stats['kept']} of {stats['read']} blocks passed")


def _window(args: argparse.Namespace, reader: MafReader) -> None:
    stats: Counter = Counter()
    _write(args, window_blocks(reader, args.length, args.step, stats), reader)
    _note(
        f"{stats['windows']} windows from {stats['read']} blocks; "
        f"{stats['skipped']} skipped without reference bases"
    )


def _merge(args: argparse.Namespace, reader: MafReader) -> None:
    options = MergeOptions(
        args.species_consensus, args.max_distance, args.max_length, args.max_gaps,
        args.min_seqs, args.reference,
    )
    stats: Counter = Counter()
    _write(args, merge_blocks(reader, options, stats), reader)
    _note(
        f"{stats['read']} blocks read, {stats['merges']} merges, {stats['written']} written, "
        f"{stats['dropped']} dropped with fewer than {args.min_seqs} sequences"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="maftools", description="Process Multiple Alignment Format (MAF) files."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", metavar="program", required=True)

    def command(name, run, help, output=("-o", "--output")):
        sub = commands.add_parser(
            name, help=help, description=help,
            formatter_class=_HelpFormatter,
        )
        sub.set_defaults(run=run)
        if output:
            sub.add_argument(*output, dest="output", default=None, metavar="FILE",
                             help="write here instead of to standard output")
        return sub

    command("describe", _describe, "Summarise a MAF file.", output=())
    command("blockdescribe", _blockdescribe, "Tabulate statistics for every block (TSV).")
    command("toBed", _to_bed, "Write the reference interval of every block as BED.",
            output=("-o", "--output", "--output_bed"))

    sub = command("toFasta", _to_fasta, "Write aligned sequences as FASTA.")
    style = sub.add_mutually_exclusive_group()
    style.add_argument("-c", "--concatenate", action="store_true",
                       help="join all blocks into one aligned record per species")
    style.add_argument("-u", "--ungapped", action="store_true", help="remove gaps")

    command("index", _index, "Index an uncompressed MAF file to speed up select.", output=())

    sub = command("select", _select, "Extract regions given in a BED file.")
    sub.add_argument("-b", "--bed", required=True, metavar="FILE", help="regions to extract")
    sub.add_argument("-l", "--left-flank", type=_at_least(0), default=0, metavar="N",
                     help="extend regions to the left by this many nucleotides")
    sub.add_argument("-r", "--right-flank", type=_at_least(0), default=0, metavar="N",
                     help="extend regions to the right by this many nucleotides")

    defaults = FilterCriteria()
    sub = command("filter", _filter, "Keep blocks that meet all criteria.")
    sub.add_argument("-g", "--max-gaps", type=_fraction, default=defaults.max_gaps, metavar="F",
                     help="maximum average gap fraction")
    sub.add_argument("-a", "--max-similarity", type=_fraction, default=defaults.max_similarity,
                     metavar="F", help="maximum identity of any sequence to the reference")
    sub.add_argument("-m", "--min-similarity", type=_fraction, default=defaults.min_similarity,
                     metavar="F", help="minimum identity of every sequence to the reference")
    sub.add_argument("-l", "--min-length", type=_at_least(0), default=defaults.min_length,
                     metavar="N", help="minimum number of columns")
    sub.add_argument("-n", "--min-seqs", type=_at_least(1), default=defaults.min_seqs,
                     metavar="N", help="minimum number of sequences")

    sub = command("species", _species, "Keep, drop or reorder species.")
    sub.add_argument("-s", "--species", action="append", default=[], metavar="NAMES",
                     help="comma-separated species or full sequence names; can be repeated")
    sub.add_argument("-f", "--list", metavar="FILE", help="file with one name per line")
    sub.add_argument("-x", "--exclude", action="store_true",
                     help="drop these species instead of keeping them")
    sub.add_argument("-a", "--require-all", action="store_true",
                     help="keep only blocks that contain all of these species")
    sub.add_argument("-r", "--reorder", action="store_true",
                     help="order sequences as given and make the first species the reference, "
                     "on the plus strand; blocks without it are dropped")
    sub.add_argument("-n", "--min-seqs", type=_at_least(1), default=1, metavar="N",
                     help="minimum number of sequences to write a block")

    sub = command("window", _window, "Cut blocks into overlapping windows.")
    sub.add_argument("-l", "--length", type=_at_least(1), default=120, metavar="N",
                     help="window length in columns; shorter blocks are kept as they are")
    sub.add_argument("-s", "--step", type=_at_least(1), default=40, metavar="N",
                     help="distance between window starts")

    merge = MergeOptions()
    sub = command("merge", _merge, "Merge consecutive blocks.")
    sub.add_argument("-s", "--species-consensus", type=_fraction, metavar="F",
                     default=merge.species_consensus,
                     help="minimum overlap of the two sets of sequence names")
    sub.add_argument("-d", "--max-distance", type=_at_least(0), default=merge.max_distance,
                     metavar="N", help="maximum number of unaligned bases between blocks")
    sub.add_argument("-l", "--max-length", type=_at_least(1), default=merge.max_length,
                     metavar="N", help="maximum number of columns of a merged block")
    sub.add_argument("-g", "--max-gaps", type=_fraction, default=merge.max_gaps, metavar="F",
                     help="drop sequences with a larger gap fraction")
    sub.add_argument("-m", "--min-seqs", type=_at_least(1), default=merge.min_seqs, metavar="N",
                     help="minimum number of sequences to merge, and to write a block")
    sub.add_argument("-r", "--no-reference", dest="reference", action="store_false",
                     help="apply --max-gaps to the first sequence as well")

    for sub in commands.choices.values():
        sub.add_argument("maf", metavar="MAF", help="input file, may be gzipped; - for standard input")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "window" and args.step > args.length:
        parser.error("--step must not exceed --length")
    if args.command == "species" and args.exclude and (args.reorder or args.require_all):
        parser.error("--exclude cannot be combined with --reorder or --require-all")
    side_inputs = (getattr(args, "bed", None), getattr(args, "list", None))
    if args.maf == STDIO and STDIO in side_inputs:
        parser.error("only one input can be read from standard input")
    try:
        with open_input(args.maf) as handle:
            args.run(args, MafReader(handle, args.maf))
        sys.stdout.flush()
    except BrokenPipeError:  # e.g. piped into head
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
    except UnicodeDecodeError:
        _note(f"error: {args.maf} or another input is not a text file")
        return 1
    except (EOFError, zlib.error) as error:
        _note(f"error: damaged gzip input: {error}")
        return 1
    except (MafError, OSError, sqlite3.Error) as error:
        _note(f"error: {error}")
        return 1
    return 0
