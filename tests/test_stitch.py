import os
from collections import Counter
from pathlib import Path

import pytest

from maftools.bed import Interval
from maftools.blocks import Block, reverse_complement
from maftools.cli import main
from maftools.mafio import format_block
from maftools.stitch import feature_name, file_names, stitch

ROOT = Path(__file__).parent.parent
MAF, BED = str(ROOT / "Example/test.maf"), str(ROOT / "Example/test.bed")


def one(blocks, start, end, strand="+", stats=None):
    (feature,) = stitch(blocks, [Interval("chr1", start, end, "f", strand)], stats)
    widths = {len(text) for text in feature.texts.values()}
    assert len(widths) == 1
    return feature


def tiled(oracle, edges, gaps=0.3, strands="+-"):
    """Blocks whose references cover the pieces between edges."""
    blocks = []
    for a, b in zip(edges, edges[1:]):
        strand = oracle.rng.choice(strands)
        start = a if strand == "+" else 5000 - b
        mask = "".join("x" + "-" * (oracle.rng.random() < gaps) for _ in range(b - a))
        rows = [oracle.row(0, start, mask, strand)]
        for i in oracle.rng.sample([1, 2, 3], 2):
            other = "".join("-" if oracle.rng.random() < gaps else "x" for _ in mask)
            rows.append(oracle.row(i, oracle.rng.randrange(4000), other,
                                   oracle.rng.choice("+-")))
        blocks.append(Block(rows))
    return blocks


@pytest.mark.parametrize("strand", "+-.")
def test_reference_is_exactly_the_feature(oracle, strand):
    truth = oracle.genomes["sp0.chr1"]
    for _ in range(10):
        edges = sorted(oracle.rng.sample(range(1000, 2000), 8))
        blocks = tiled(oracle, edges)
        oracle.rng.shuffle(blocks)
        start = oracle.rng.randrange(edges[0], edges[-1] - 1)
        end = oracle.rng.randrange(start + 1, edges[-1])
        feature = one(blocks, start, end, strand)
        want = truth[start:end]
        if strand == "-":
            want = reverse_complement(want)
        assert feature.texts["sp0"].replace("-", "") == want
        assert list(feature.texts)[0] == "sp0" and feature.covered == end - start


def test_other_species_keep_their_alignment_to_the_reference(oracle):
    """Every aligned pair of bases in the result occurs in the blocks, and all of them."""
    blocks = tiled(oracle, [1000, 1080, 1200, 1230])

    def pairs(texts, position):
        out = []
        for column in zip(*texts):
            if column[0] != "-":
                out.append((position, column[1:]))
                position += 1
        return out

    names = ["sp1", "sp2", "sp3"]
    want = []
    for block in blocks:
        if block.reference.strand == "-":
            block = block.reverse_complement()
        texts = {row.species: row.text for row in block.rows}
        empty = "-" * block.columns
        want += pairs([texts["sp0"], *(texts.get(n, empty) for n in names)],
                      block.reference.start)
    feature = one(blocks, 1000, 1230)
    empty = "-" * len(feature.texts["sp0"])
    got = pairs([feature.texts["sp0"], *(feature.texts.get(n, empty) for n in names)], 1000)
    assert got == want


def test_uncovered_positions_are_n(oracle):
    blocks = [Block([oracle.row(i, s, "x" * 50) for i in (0, 1)]) for s in (1000, 1080)]
    stats = Counter()
    feature = one(blocks, 990, 1150, stats=stats)
    truth = oracle.genomes["sp0.chr1"]
    assert feature.texts["sp0"] == (
        "N" * 10 + truth[1000:1050] + "N" * 30 + truth[1080:1130] + "N" * 20)
    assert feature.texts["sp1"].count("-") == 60 and feature.covered == 100
    assert stats["partial"] == 1 and stats["matched"] == 1


def test_overlapping_blocks_are_cut(oracle):
    a = Block([oracle.row(i, 1000, "x" * 100) for i in (0, 1)])
    b = Block([oracle.row(i, 1060, "x" * 100) for i in (0, 2)])
    inside = Block([oracle.row(i, 1020, "x" * 10) for i in (0, 3)])
    stats = Counter()
    feature = one([b, inside, a], 1000, 1160, stats=stats)
    assert feature.texts["sp0"] == oracle.genomes["sp0.chr1"][1000:1160]
    assert feature.texts["sp2"] == "-" * 100 + oracle.genomes["sp2.chr1"][1100:1160]
    assert "sp3" not in feature.texts and stats["overlaps"] == 2 and feature.covered == 160


def test_identical_blocks_are_used_once(oracle):
    block = Block([oracle.row(i, 1000, "x" * 100) for i in (0, 1)])
    stats = Counter()
    feature = one([block, block], 1010, 1060, stats=stats)
    assert feature.texts["sp0"] == oracle.genomes["sp0.chr1"][1010:1060]
    assert feature.covered == 50 and stats["overlaps"] == 1


def test_features_keep_their_order_whatever_the_order_of_blocks(oracle):
    blocks = [Block([oracle.row(0, s, "x" * 50)]) for s in (2000, 1000)]
    intervals = [Interval("chr1", 1010, 1020, "a"), Interval("chr1", 2010, 2020, "b")]
    assert [f.name for f in stitch(blocks, intervals)] == ["a", "b"]


def test_first_copy_of_a_species_is_used(oracle):
    block = Block([oracle.row(0, 0, "x" * 9), oracle.row(1, 0, "x" * 9), oracle.row(1, 50, "x" * 9)])
    stats = Counter()
    assert one([block], 0, 9, stats=stats).texts["sp1"] == oracle.genomes["sp1.chr1"][:9]
    assert stats["duplicates"] == 1


def test_end_is_limited_to_the_sequence(oracle):
    block = Block([oracle.row(0, 4950, "x" * 50)])
    feature = one([block], 4960, 5100)
    assert feature.end == 5000 and len(feature.texts["sp0"]) == 40 and feature.covered == 40


def test_order_and_unmatched(oracle):
    block = Block([oracle.row(0, 1000, "x" * 100)])
    stats = Counter()
    found = list(stitch([block], [Interval("chr1", 1050, 1060, "b"), Interval("chr1", 1, 2, "no"),
                                  Interval("chr1", 1010, 1020, "a")], stats))
    assert [f.name for f in found] == ["b", "a"]
    assert (stats["intervals"], stats["matched"], stats["partial"]) == (3, 2, 0)


def test_names():
    assert feature_name(Interval("c", 1, 2, "tRNA", "-")) == "tRNA"
    assert feature_name(Interval("c", 1, 2, "my gene")) == "my_gene"
    assert feature_name(Interval("c", 1, 2, ".", "-")) == "c_1_2_-"
    names = list(file_names(["a", "A", "a", "../../etc/x", "..", "", "a b/c", "x" * 500]))
    assert names == ["a", "A_2", "a_3", "_.._etc_x", "_", "__2", "a_b_c", "x" * 200]
    assert not any(os.sep in name or name.startswith(".") for name in names)


def run(capsys, *argv):
    code = main([str(a) for a in argv])
    out, err = capsys.readouterr()
    return code, out, err


def test_cli(capsys, tmp_path):
    code, out, err = run(capsys, "stitch", "-b", BED, MAF)
    assert code == 0 and "3 of 5 intervals matched, 2 only in part" in err
    name = "NC_037638.1_1055955_1056026_+"
    assert f">{name}|NC_037638 NC_037638.1:1055955-1056026(+)\n" in out
    folder = tmp_path / "new" / "folder"
    assert run(capsys, "stitch", "-b", BED, "-d", folder, MAF)[0] == 0
    files = sorted(os.listdir(folder))
    assert len(files) == 3 and name + ".fa" in files
    assert (folder / (name + ".fa")).read_text().startswith(">NC_037638 NC_037638.1:")
    with pytest.raises(SystemExit):
        main(["stitch", "-b", BED, "-d", "x", "-o", "y", MAF])


def test_cli_with_index_and_flanks(capsys, tmp_path):
    maf = tmp_path / "x.maf"
    maf.write_text(Path(MAF).read_text())
    plain = run(capsys, "stitch", "-l", "40", "-r", "9", "-b", BED, maf)
    run(capsys, "index", maf)
    assert run(capsys, "stitch", "-l", "40", "-r", "9", "-b", BED, maf) == plain


def test_cli_flanks(capsys):
    lengths = []
    for flanks in ((), ("-l", "3"), ("-r", "2"), ("-l", "3", "-r", "2")):
        out = run(capsys, "stitch", *flanks, "-b", BED, MAF)[1]
        reference = out.split(">")[1].split("\n", 1)[1].replace("\n", "").replace("-", "")
        lengths.append(len(reference))
    assert lengths == [71, 74, 73, 76]


def test_flanks_stay_within_the_sequence(capsys, tmp_path, oracle):
    maf, bed = tmp_path / "x.maf", tmp_path / "x.bed"
    maf.write_text(format_block(Block([oracle.row(0, 0, "x" * 50)])))
    bed.write_text("chr1\t2\t10\n")
    out = run(capsys, "stitch", "-l", "5", "-b", bed, maf)[1]
    assert "sp0.chr1:0-10(.)" in out and out.split("\n")[1] == oracle.genomes["sp0.chr1"][:10]


def test_cli_select_stranded(capsys, tmp_path, oracle):
    maf, bed = tmp_path / "x.maf", tmp_path / "x.bed"
    block = Block([oracle.row(0, 100, "x" * 50), oracle.row(1, 100, "x" * 50)])
    maf.write_text(format_block(block))
    bed.write_text("chr1\t110\t120\tx\t0\t-\n")
    assert " + " in run(capsys, "select", "-b", bed, maf)[1]
    out = run(capsys, "select", "-s", "-b", bed, maf)[1]
    assert " - " in out and " + " not in out


def test_outputs_do_not_overwrite_inputs(capsys, tmp_path):
    maf, bed = tmp_path / "f.fa", tmp_path / "x.bed"
    maf.write_text(Path(MAF).read_text())
    bed.write_text("NC_037638.1\t1055955\t1056026\tf\t0\t+\n")
    assert run(capsys, "stitch", "-b", bed, "-d", tmp_path, maf)[0] == 1
    assert maf.read_text() == Path(MAF).read_text()
