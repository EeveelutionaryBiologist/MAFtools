"""Command line interface."""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from collections.abc import Iterable, Iterator, Sequence
from contextlib import contextmanager
from typing import TextIO

from . import __version__
from .bed import read_bed
from .blocks import Block
from .mafio import MafError, MafReader, write_maf
from .operations import (
    FilterCriteria, MergeOptions, filter_blocks, merge_blocks, select_blocks, window_blocks,
)
from .reports import bed_lines, block_table, describe


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
    if args.output is None:
        yield sys.stdout
        return
    for path in (args.maf, getattr(args, "bed", None)):
        if path and os.path.exists(args.output) and os.path.samefile(path, args.output):
            raise MafError(f"{args.output} is also an input file")
    with open(args.output, "w") as handle:
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


def _select(args: argparse.Namespace, reader: MafReader) -> None:
    with open(args.bed) as handle:
        intervals = read_bed(handle, args.bed)
    stats: Counter = Counter()
    blocks = select_blocks(reader, intervals, args.left_flank, args.right_flank, stats)
    written = _write(args, blocks, reader)
    _note(f"{written} blocks written; {stats['matched']} of {stats['intervals']} intervals matched")


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
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        sub.add_argument("maf", metavar="MAF", help="input file")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "window" and args.step > args.length:
        parser.error("--step must not exceed --length")
    try:
        with open(args.maf) as handle:
            args.run(args, MafReader(handle, args.maf))
        sys.stdout.flush()
    except BrokenPipeError:  # e.g. piped into head
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
    except UnicodeDecodeError:
        _note(f"error: {args.maf} or another input is not a text file")
        return 1
    except (MafError, OSError) as error:
        _note(f"error: {error}")
        return 1
    return 0
