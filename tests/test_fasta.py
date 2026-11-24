from collections import Counter

from maftools.blocks import Block, Row
from maftools.reports import fasta_blocks, fasta_concatenated


def parse(lines):
    out = {}
    for line in lines:
        if line.startswith(">"):
            name = line[1:]
            assert name not in out
            out[name] = ""
        else:
            assert len(line) <= 60
            out[name] += line
    return out


def test_blocks(oracle):
    block = oracle.random_block(columns=130)
    records = parse(fasta_blocks([block, block]))
    assert len(records) == 8
    for row, (name, text) in zip(block.rows, records.items()):
        start, end = row.plus_interval()
        assert name == f"{row.src}:{start}-{end}({row.strand}) block=1" and text == row.text


def test_ungapped_records_are_the_named_interval(oracle):
    block = oracle.random_block(columns=130)
    complement = str.maketrans("ACGT", "TGCA")
    for row, (name, text) in zip(block.rows, parse(fasta_blocks([block], True)).items()):
        start, end = row.plus_interval()
        truth = oracle.genomes[row.src][start:end]
        assert text == (truth if row.strand == "+" else truth.translate(complement)[::-1])


def row(src, text):
    return Row(src, 0, len(text.replace("-", "")), "+", 100, text)


def test_concatenated_pads_missing_species():
    blocks = [
        Block([row("a.1", "AC"), row("b.1", "GT")]),
        Block([row("a.2", "TTT"), row("c.1", "CCC")]),
        Block([row("b.9", "G")]),
    ]
    assert parse(fasta_concatenated(blocks)) == {
        "a": "ACTTT-", "b": "GT---G", "c": "--CCC-",
    }


def test_concatenated_uses_first_copy():
    stats = Counter()
    block = Block([row("a.1", "AC"), row("b.1", "GT"), row("b.2", "TT")])
    assert parse(fasta_concatenated([block], stats)) == {"a": "AC", "b": "GT"}
    assert stats["duplicates"] == 1


def test_nothing():
    assert list(fasta_blocks([])) == list(fasta_concatenated([])) == []
