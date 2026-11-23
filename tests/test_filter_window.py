from collections import Counter

import pytest

from maftools.blocks import Block
from maftools.metrics import pairwise_identity
from maftools.operations import FilterCriteria, filter_blocks, window_blocks, window_starts

LOOSE = dict(max_gaps=1, min_similarity=0, max_similarity=1, min_length=0, min_seqs=1)


def accepts(block, **kw):
    return FilterCriteria(**{**LOOSE, **kw}).accepts(block)


def test_min_length_is_a_minimum(oracle):
    block = oracle.random_block(columns=100)
    assert [accepts(block, min_length=n) for n in (40, 99, 100, 101)] == [True, True, True, False]


def test_min_seqs_and_gaps(oracle):
    rows = [oracle.row(0, 0, "x" * 10), oracle.row(1, 0, "x" * 5 + "-" * 5)]
    block = Block(rows)  # mean gap fraction 0.25
    assert accepts(block, min_seqs=2) and not accepts(block, min_seqs=3)
    assert accepts(block, max_gaps=0.25) and not accepts(block, max_gaps=0.24)


def test_similarity_applies_to_every_row(oracle):
    a, b, c = (oracle.row(i, 0, "x" * 10) for i in range(3))
    b.text = a.text
    c.text = a.text[:5] + a.text[5:].translate(str.maketrans("ACGT", "CATG"))
    block = Block([a, b, c])
    assert accepts(block, min_similarity=0.5) and not accepts(block, min_similarity=0.6)
    assert not accepts(block, max_similarity=0.9)


def test_filter_passes_blocks_through_untouched(oracle):
    block = oracle.random_block()
    block.empties.append("rn.chr3 5 10 + 300 I")
    stats = Counter()
    assert list(filter_blocks([block], FilterCriteria(**LOOSE), stats))[0] is block
    assert stats == {"read": 1, "kept": 1}


def test_identity():
    assert pairwise_identity("AC-G", "ac-T") == pytest.approx(2 / 3)
    assert pairwise_identity("--", "--") == 0
    with pytest.raises(ValueError):
        pairwise_identity("A", "AC")


def test_windows_are_true_to_genomes(oracle):
    block = oracle.random_block(columns=250, gaps=0.5)
    windows = list(window_blocks([block], 40, 15))
    assert all(w.columns == 40 for w in windows)
    for window in windows:
        oracle.check(window)


def test_short_block_is_kept_as_it_is(oracle):
    block = oracle.random_block(columns=40)
    assert list(window_blocks([block], 40, 20)) == [block]


def test_window_starts_cover_the_block():
    assert window_starts(100, 40, 20) == [0, 20, 40, 60]
    assert window_starts(105, 40, 20) == [0, 20, 40, 60, 65]
    for n in range(1, 60):
        covered = {c for s in window_starts(n, 10, 7) for c in range(s, min(n, s + 10))}
        assert covered == set(range(n))


def test_window_without_reference_bases_is_skipped(oracle):
    block = Block([oracle.row(0, 0, "x" * 10 + "-" * 20), oracle.row(1, 0, "x" * 30)])
    stats = Counter()
    assert len(list(window_blocks([block], 10, 10, stats))) == 1
    assert stats["skipped"] == 2


def test_window_step_is_validated(oracle):
    for step in (0, 11):
        with pytest.raises(ValueError):
            list(window_blocks([], 10, step))
