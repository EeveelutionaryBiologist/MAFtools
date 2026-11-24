from collections import Counter

import pytest

from maftools.blocks import Block
from maftools.operations import dedup_blocks


def dedup(block, keep="best", min_seqs=1, stats=None):
    return list(dedup_blocks([block], keep, min_seqs, stats))


def starts(block):
    return [(row.species, row.start) for row in block.rows]


@pytest.fixture
def block(oracle):
    """sp1 has three copies; the second one is identical to the reference."""
    ref = oracle.genomes["sp0.chr1"][100:160]
    genome = oracle.genomes["sp1.chr1"]
    oracle.genomes["sp1.chr1"] = genome[:2000] + ref + genome[2060:]
    rows = [(0, 100), (1, 500), (2, 100), (1, 2000), (1, 3000)]
    return Block([oracle.row(i, s, "x" * 60) for i, s in rows], {"score": "5"})


def test_best_keeps_the_copy_most_similar_to_the_reference(oracle, block):
    (out,) = dedup(block)
    assert starts(out) == [("sp0", 100), ("sp2", 100), ("sp1", 2000)]
    oracle.check(out)
    assert "score" not in out.attrs


def test_first(oracle, block):
    (out,) = dedup(block, "first")
    assert starts(out) == [("sp0", 100), ("sp1", 500), ("sp2", 100)]
    oracle.check(out)


def test_none(oracle, block):
    (out,) = dedup(block, "none")
    assert starts(out) == [("sp0", 100), ("sp2", 100)]
    oracle.check(out)


def test_ties_go_to_the_first_copy(oracle):
    a, b = oracle.row(1, 500, "x" * 20), oracle.row(1, 900, "x" * 20)
    a.text = b.text = "A" * 20
    (out,) = dedup(Block([oracle.row(0, 100, "x" * 20), a, b]))
    assert starts(out) == [("sp0", 100), ("sp1", 500)]


def test_gaps_count_against_a_copy(oracle):
    ref = oracle.genomes["sp0.chr1"][100:140]
    for start in (500, 900):
        genome = oracle.genomes["sp1.chr1"]
        oracle.genomes["sp1.chr1"] = genome[:start] + ref + genome[start + 40:]
    fragment, whole = oracle.row(1, 500, "x" * 10 + "-" * 30), oracle.row(1, 900, "x" * 40)
    (out,) = dedup(Block([oracle.row(0, 100, "x" * 40), fragment, whole]))
    assert starts(out) == [("sp0", 100), ("sp1", 900)]


def test_copies_on_other_sequences_of_a_species_count(oracle):
    rows = [oracle.row(0, 0, "x" * 20), oracle.row(1, 0, "x" * 20), oracle.row(1, 50, "x" * 20)]
    rows[2].src = "sp1.chr7"
    (out,) = dedup(Block(rows), "first")
    assert [row.src for row in out.rows] == ["sp0.chr1", "sp1.chr1"]


@pytest.mark.parametrize("keep", ["best", "first"])
def test_reference_is_kept_over_its_copies(oracle, keep):
    rows = [oracle.row(0, 100, "x" * 20), oracle.row(0, 100, "x" * 20), oracle.row(1, 0, "x" * 20)]
    rows[0].text = "N" * 20
    (out,) = dedup(Block(rows), keep)
    assert out.reference.text == "N" * 20 and len(out.rows) == 2


def test_none_drops_blocks_with_a_duplicated_reference(oracle):
    rows = [oracle.row(0, 100, "x" * 20), oracle.row(1, 0, "x" * 20), oracle.row(0, 900, "x" * 20)]
    stats = Counter()
    assert dedup(Block(rows), "none", stats=stats) == []
    assert stats == {"read": 1, "affected": 1, "dropped": 1}


def test_gap_columns_are_removed(oracle):
    rows = [oracle.row(0, 0, "xx--xx"), oracle.row(1, 0, "xx--xx"), oracle.row(1, 90, "--xx--")]
    (out,) = dedup(Block(rows), "first")
    assert out.columns == 4
    oracle.check(out)


def test_block_without_duplicates_is_unchanged(oracle):
    block = oracle.random_block()
    for keep in ("best", "first", "none"):
        (out,) = dedup(block, keep)
        assert out is block


def test_random_blocks(oracle):
    for _ in range(20):
        rows = [row for _ in range(3) for row in oracle.random_block(gaps=0.4).rows]
        oracle.rng.shuffle(rows)
        for keep in ("best", "first"):
            (out,) = dedup(Block(rows), keep)
            oracle.check(out)
            assert sorted(row.species for row in out.rows) == ["sp0", "sp1", "sp2", "sp3"]
            assert out.reference.start == rows[0].start
            assert [x for x in dedup(out, keep)] == [out]


def test_min_seqs_and_stats(oracle, block):
    stats = Counter()
    assert len(dedup(block, min_seqs=3, stats=stats)) == 1
    assert dedup(block, min_seqs=4, stats=stats) == []
    assert stats == {"read": 2, "affected": 2, "removed": 2, "kept": 1, "dropped": 1}


def test_unknown_mode():
    with pytest.raises(ValueError):
        list(dedup_blocks([], "worst"))
