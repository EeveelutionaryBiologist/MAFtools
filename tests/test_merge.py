from collections import Counter

import pytest

from maftools.blocks import Block
from maftools.operations import MergeOptions, merge_blocks


def merge(*blocks, **kw):
    return list(merge_blocks(blocks, MergeOptions(**{"max_gaps": 1.0, **kw})))


def full(oracle, starts, columns=100, species=None):
    species = range(len(starts)) if species is None else species
    return Block([oracle.row(i, starts[i], "x" * columns) for i in species])


def test_split_blocks_are_restored(oracle):
    whole = oracle.random_block(columns=300, gaps=0.3).without_gap_columns()
    cuts = [whole.sliced(a, b) for a, b in ((0, 90), (90, 200), (200, whole.columns))]
    (restored,) = merge(*cuts)
    assert restored.rows == whole.rows
    oracle.check(restored)


def test_rows_are_paired_by_name_not_by_order(oracle):
    a = full(oracle, [1000] * 4)
    b = full(oracle, [1100] * 4, species=[0, 3, 2, 1])
    (merged,) = merge(a, b)
    oracle.check(merged)
    assert [r.src for r in merged.rows] == [f"sp{i}.chr1" for i in range(4)]


def test_unequal_sizes_and_distances(oracle):
    a = Block([oracle.row(i, 1000, "x" * (100 - 5 * i) + "-" * 5 * i) for i in range(4)])
    b = full(oracle, [1100, 1099, 1110, 1085])
    (merged,) = merge(a, b, max_distance=50)
    oracle.check(merged, bridged=True)
    assert [r.size for r in merged.rows] == [200, 199, 210, 185]
    assert [r.text.count("N") for r in merged.rows] == [0, 4, 20, 0]


def test_inputs_are_not_modified(oracle):
    a, b = full(oracle, [1000] * 3), full(oracle, [1100] * 3)
    before = repr((a, b))
    merge(a, b)
    assert repr((a, b)) == before


@pytest.mark.parametrize(
    "starts, why",
    [
        ([1100, 2000, 2000], "only the reference continues"),
        ([500, 500, 500], "out of order"),
        ([1050, 1050, 1050], "overlapping"),
        ([1101, 1101, 1101], "too far apart"),
    ],
)
def test_unmergeable_blocks_are_both_kept(oracle, starts, why):
    a, b = full(oracle, [1000] * 3), full(oracle, starts)
    assert merge(a, b) == [a, b], why


def test_strand_and_reference_must_agree(oracle):
    a = full(oracle, [1000] * 3)
    flipped = Block([oracle.row(i, 1100, "x" * 50, "-" if i == 0 else "+") for i in range(3)])
    other = full(oracle, [1100] * 3, species=[1, 0, 2])
    assert len(merge(a, flipped)) == 2 and len(merge(a, other)) == 2


def test_row_on_other_strand_is_not_joined(oracle):
    a = full(oracle, [1000] * 3)
    b = Block([oracle.row(i, 1100, "x" * 50, "-" if i == 2 else "+") for i in range(3)])
    (merged,) = merge(a, b)
    assert [r.src for r in merged.rows] == ["sp0.chr1", "sp1.chr1"]


def test_max_length_includes_the_bridge(oracle):
    a, b = full(oracle, [1000] * 3), full(oracle, [1110] * 3)
    assert len(merge(a, b, max_distance=10, max_length=209)) == 2
    assert len(merge(a, b, max_distance=10, max_length=210)) == 1


def test_consensus_and_min_seqs(oracle):
    a, b = full(oracle, [1000] * 4), full(oracle, [1100] * 4, species=[0, 1])
    assert len(merge(a, b, species_consensus=0.67)) == 2  # 2 * 2 / (4 + 2)
    assert len(merge(a, b, species_consensus=0.66)) == 1
    two = full(oracle, [1000] * 2), full(oracle, [1100] * 2)
    assert len(merge(*two, min_seqs=2)) == 1
    assert merge(*two, min_seqs=3) == []  # neither merged nor large enough to keep


def test_blocks_sharing_too_few_rows_are_kept_apart(oracle):
    a = full(oracle, [1000] * 4, species=[0, 1, 2])
    b = full(oracle, [1100] * 4, species=[0, 1, 3])
    assert merge(a, b, species_consensus=0.5, min_seqs=3) == [a, b]


def test_last_block_below_min_seqs(oracle):
    a, b = full(oracle, [1000] * 3), full(oracle, [3000], species=[0])
    assert merge(a, b) == [a]


def test_duplicated_source_uses_the_nearest_copy(oracle):
    a = Block([oracle.row(0, 1000, "x" * 50), oracle.row(1, 1000, "x" * 50),
               oracle.row(1, 2000, "x" * 50)])
    b = Block([oracle.row(0, 1050, "x" * 50), oracle.row(1, 2050, "x" * 50),
               oracle.row(1, 1050, "x" * 50)])
    (merged,) = merge(a, b)
    oracle.check(merged)
    assert [r.start for r in merged.rows] == [1000, 1000, 2000]


def test_duplicated_source_is_paired_once_and_nearest_first(oracle):
    far, near = oracle.row(1, 990, "x" * 50), oracle.row(1, 1000, "x" * 50)
    a = Block([oracle.row(0, 1000, "x" * 50), far, near])
    b = Block([oracle.row(i, 1050, "x" * 50) for i in (0, 1)])
    (merged,) = merge(a, b, max_distance=10)
    oracle.check(merged, bridged=True)
    assert [r.start for r in merged.rows] == [1000, 1000] and len(merged.rows) == 2


def test_references_are_paired_even_if_another_copy_is_closer(oracle):
    a = Block([oracle.row(i, 1000, "x" * 50) for i in (0, 1)])
    b = Block([oracle.row(0, 1055, "x" * 50), oracle.row(0, 1052, "x" * 50),
               oracle.row(1, 1050, "x" * 50)])
    (merged,) = merge(a, b, max_distance=10)
    oracle.check(merged, bridged=True)
    assert merged.reference.size == 105


def test_gappy_rows_and_resulting_gap_columns_are_removed(oracle):
    rows = [oracle.row(0, 0, "xx--xxxx"), oracle.row(1, 0, "xx--xxxx"),
            oracle.row(2, 0, "--xx----")]
    (out,) = merge(Block(rows), max_gaps=0.5)
    assert out.columns == 6 and len(out.rows) == 2
    oracle.check(out)


def test_reference_is_exempt_from_max_gaps(oracle):
    rows = [oracle.row(0, 0, "x-------"), *(oracle.row(i, 0, "x" * 8) for i in (1, 2))]
    kept, dropped = merge(Block(rows), max_gaps=0.5), merge(
        Block(rows), max_gaps=0.5, keep_reference=False)
    assert kept[0].reference.src == "sp0.chr1" and dropped[0].reference.src == "sp1.chr1"


def test_stats(oracle):
    stats = Counter()
    blocks = [full(oracle, [1000] * 3), full(oracle, [1100] * 3), full(oracle, [3000], species=[0])]
    list(merge_blocks(blocks, MergeOptions(), stats))
    assert stats == {"read": 3, "merges": 1, "written": 1, "dropped": 1}


def test_no_blocks():
    assert merge() == []
