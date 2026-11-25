import pytest
from io import StringIO

from maftools.bed import Interval, IntervalIndex, chrom_names, read_bed
from maftools.blocks import Block
from maftools.mafio import MafError
from maftools.operations import select_blocks


def select(block, start, end, chrom="sp0.chr1", **kw):
    return list(select_blocks([block], [Interval(chrom, start, end)], **kw))


def test_random_regions_are_true_to_genomes(oracle):
    for _ in range(40):
        block = oracle.random_block()
        lo, hi = block.reference.plus_interval()
        start = oracle.rng.randrange(lo, hi)
        end = oracle.rng.randrange(start + 1, hi + 1)
        (part,) = select(block, start, end)
        oracle.check(part)
        assert part.reference.plus_interval() == (start, end)


def test_whole_block_is_unchanged(oracle):
    block = oracle.random_block(gaps=0.4)
    assert select(block, 0, 5000) == [block]


def test_adjacent_intervals_partition_columns(oracle):
    block = oracle.random_block(gaps=0.4, strands="+")
    lo, hi = block.reference.plus_interval()
    mid = (lo + hi) // 2
    (left,), (right,) = select(block, lo, mid), select(block, mid, hi)
    assert left.columns + right.columns == block.columns


def test_single_base(oracle):
    block = oracle.random_block()
    lo, _ = block.reference.plus_interval()
    (part,) = select(block, lo + 10, lo + 11)
    assert part.reference.size == 1


def test_gap_only_row_is_dropped_but_block_kept(oracle):
    rows = [oracle.row(0, 1000, "x" * 100), oracle.row(1, 1000, "x" * 10 + "-" * 90)]
    (part,) = select(Block(rows), 1050, 1060)
    assert [r.src for r in part.rows] == ["sp0.chr1"]


def test_flanks(oracle):
    block = Block([oracle.row(0, 1000, "x" * 100), oracle.row(1, 1000, "x" * 100)])
    (part,) = select(block, 1050, 1060, left=5, right=7)
    assert part.reference.plus_interval() == (1045, 1067)


@pytest.mark.parametrize("strand", "+-.")
def test_stranded(oracle, strand):
    for _ in range(10):
        block = oracle.random_block()
        lo, hi = block.reference.plus_interval()
        interval = Interval("chr1", lo + 10, hi - 20, "x", strand)
        (part,) = select_blocks([block], [interval], stranded=True)
        (plain,) = select_blocks([block], [interval])
        oracle.check(part)
        assert part.reference.strand == (block.reference.strand if strand == "." else strand)
        assert [r.plus_interval() for r in part.rows] == [r.plus_interval() for r in plain.rows]


def test_no_overlap(oracle):
    block = Block([oracle.row(0, 1000, "x" * 100)])
    assert select(block, 900, 1000) == [] and select(block, 1100, 1200) == []


def test_slicing_drops_stale_score(oracle):
    block = oracle.random_block(strands="+")
    (part,) = select(block, 1010, 1020)
    assert "score" not in part.attrs and block.attrs["score"] == "1234.0"


@pytest.mark.parametrize(
    "src, bed, matches",
    [
        ("hg38.chr1", "hg38.chr1", True),
        ("hg38.chr1", "chr1", True),
        ("hg38.chr10", "chr1", False),
        ("hg38.chr10", "hg38.chr1", False),
        ("NC_037638.1", "NC_037638.1", True),
        ("NC_037638.1", "1", False),
        ("chr1", "chr1", True),
    ],
)
def test_chromosome_names(src, bed, matches):
    assert (bed in chrom_names(src)) is matches


def test_index_agrees_with_brute_force(oracle):
    rng = oracle.rng
    ivs = [Interval("c", s, s + rng.randrange(0, 300)) for s in rng.choices(range(2000), k=200)]
    index = IntervalIndex(ivs)
    for _ in range(50):
        lo = rng.randrange(2000)
        hi = lo + rng.randrange(1, 200)
        want = [(i, s, e) for i, (_, s, e, *_) in enumerate(ivs) if max(s, lo) < min(e, hi)]
        assert index.overlapping("c", lo, hi) == want


def test_bed_parsing():
    text = "# c\ntrack name=x\nbrowser position\n\nchr1\t5\t9\tname\t0\t-\n"
    assert read_bed(StringIO(text)) == [Interval("chr1", 5, 9, "name", "-")]
    assert read_bed(StringIO("chr1\t5\t9\n")) == [Interval("chr1", 5, 9, "", ".")]
    assert read_bed(StringIO("chr1 5 9 n 0 +\n")) == [Interval("chr1", 5, 9, "n", "+")]
    assert read_bed(StringIO("chr1\t5\t9\t.\t-\n")) == [Interval("chr1", 5, 9, ".", "-")]
    assert read_bed(StringIO("chr1\t5\t9\t.\t100\n"), stranded=True)[0].strand == "."
    assert read_bed(StringIO("chr1\t5\t9\tn\t-\t+\n"))[0].strand == "+"
    spaced = "track name=a b\tc\nchr1\t5\t9\tmy gene\t0\t-\n\t\n"
    assert read_bed(StringIO(spaced), stranded=True) == [Interval("chr1", 5, 9, "my gene", "-")]
    for odd in ("1", "", "x"):
        line = f"chr1\t5\t9\tn\t0\t{odd}\tmore\n"
        assert read_bed(StringIO(line))[0].strand == "."
        with pytest.raises(MafError, match="<bed>:1: invalid strand"):
            read_bed(StringIO(line), stranded=True)
    for bad in ("chr1\t5\n", "chr1\tx\t9\n", "chr1\t9\t5\n", "chr1\t-1\t5\n"):
        with pytest.raises(MafError, match="<bed>:1"):
            read_bed(StringIO(bad))
