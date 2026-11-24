import os

import pytest

from maftools.bed import Interval
from maftools.blocks import Block
from maftools.index import MafIndex, build_index, index_path
from maftools.mafio import MafError, MafReader, format_block
from maftools.operations import select_blocks


def write(path, blocks, newline="\n"):
    text = "##maf version=1\n# comment\n\n"
    for i, block in enumerate(blocks):
        text += format_block(block) + ("# between\n\n" if i % 3 == 0 else "")
        text += "a score=1\n\n" if i % 5 == 0 else ""  # block without sequences
    with open(path, "w", newline="") as handle:
        handle.write(text.replace("\n", newline).rstrip())


def random_blocks(oracle, n=60):
    out = []
    for _ in range(n):
        block = oracle.random_block(start=oracle.rng.randrange(0, 4000),
                                    columns=oracle.rng.randrange(5, 300))
        if oracle.rng.random() < 0.3:
            block.rows[0].src = "sp0.chr10"
        out.append(block)
    return out


def open_blocks(path):
    with open(path) as handle:
        yield from MafReader(handle)


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_indexed_select_equals_full_scan(oracle, tmp_path, newline):
    maf = str(tmp_path / "x.maf")
    blocks = random_blocks(oracle)
    write(maf, blocks, newline)
    assert build_index(maf) == len(blocks)
    index = MafIndex(maf)
    assert index.is_current()
    for n in (1, 3, 40):
        intervals = []
        for _ in range(n):
            start = oracle.rng.randrange(0, 4900)
            name = oracle.rng.choice(["chr1", "sp0.chr1", "chr10", "chr2"])
            intervals.append(Interval(name, start, start + oracle.rng.randrange(0, 300)))
        scan = list(select_blocks(open_blocks(maf), intervals))
        assert list(select_blocks(index.blocks(intervals), intervals)) == scan
    index.close()


def test_only_overlapping_blocks_are_read(oracle, tmp_path):
    maf = str(tmp_path / "x.maf")
    blocks = [Block([oracle.row(0, s, "x" * 100)]) for s in range(0, 2000, 100)]
    write(maf, blocks)
    build_index(maf)
    index = MafIndex(maf)
    assert len(index.offsets([Interval("chr1", 500, 600)])) == 1
    assert len(index.offsets([Interval("chr1", 599, 601)])) == 2
    assert index.offsets([Interval("chr2", 0, 5000)]) == []
    index.close()


def test_touching_blocks_are_not_read(oracle, tmp_path):
    maf = str(tmp_path / "x.maf")
    write(maf, [Block([oracle.row(0, s, "x" * n)]) for s, n in ((0, 900), (900, 100))])
    build_index(maf)
    index = MafIndex(maf)
    assert len(index.offsets([Interval("chr1", 1000, 1100)])) == 0
    assert len(index.offsets([Interval("chr1", 900, 1100)])) == 1
    index.close()


def test_changed_file_is_noticed(oracle, tmp_path):
    maf = str(tmp_path / "x.maf")
    write(maf, random_blocks(oracle, 5))
    build_index(maf)
    info = os.stat(maf)
    os.utime(maf, ns=(info.st_atime_ns, info.st_mtime_ns + 1))
    index = MafIndex(maf)
    assert not index.is_current()
    index.close()
    build_index(maf)
    index = MafIndex(maf)
    assert index.is_current()
    index.close()


def test_failed_build_leaves_nothing_behind(tmp_path):
    maf = tmp_path / "x.maf"
    maf.write_text("a\ns x 0 9 + 9 AC\n")
    with pytest.raises(MafError):
        build_index(str(maf))
    assert os.listdir(tmp_path) == ["x.maf"]


def test_empty_file(tmp_path):
    maf = tmp_path / "x.maf"
    maf.write_text("##maf version=1\n")
    assert build_index(str(maf)) == 0
    index = MafIndex(str(maf))
    assert index.is_current() and index.offsets([Interval("chr1", 0, 9)]) == []
    index.close()


def test_gzip_is_refused(tmp_path):
    maf = tmp_path / "x.maf"
    maf.write_bytes(b"\x1f\x8b rest")
    with pytest.raises(MafError, match="uncompressed"):
        build_index(str(maf))
    assert not os.path.exists(index_path(str(maf)))
