import os
from collections import Counter
from io import StringIO

import pytest

from maftools.blocks import Block, Row, reverse_complement
from maftools.cli import main
from maftools.fasta import Fasta, Genomes
from maftools.mafio import MafError, MafReader, format_block
from maftools.operations import MergeOptions, merge_blocks, window_blocks
from maftools.verify import verify


def write_fasta(path, sequences, width=60, newline="\n"):
    with open(path, "w", newline="") as handle:
        for name, text in sequences.items():
            handle.write(f">{name} description{newline}")
            for i in range(0, len(text), width):
                handle.write(text[i : i + width] + newline)
    return str(path)


@pytest.fixture
def genomes(oracle, tmp_path):
    for src, text in oracle.genomes.items():
        write_fasta(tmp_path / (src.split(".")[0] + ".fa"), {"chr1": text})
    found = Genomes()
    assert found.add_directory(str(tmp_path)) == 5
    yield found
    found.close()


def problems(text, genomes):
    stats, unknown = Counter(), Counter()
    if isinstance(text, Block):
        text = format_block(text)
    out = list(verify(MafReader(StringIO(text)), genomes, stats, unknown))
    assert stats["problems"] == len(out)
    return out, stats, unknown


def kinds(block, genomes):
    return [p.problem for p in problems(block, genomes)[0]]


# --- FASTA access -----------------------------------------------------------

@pytest.mark.parametrize("width, newline", [(60, "\n"), (7, "\r\n"), (100000, "\n")])
def test_fetch(oracle, tmp_path, width, newline):
    path = write_fasta(tmp_path / "x.fa", oracle.genomes, width, newline)
    fasta = Fasta(path)
    for src, text in oracle.genomes.items():
        assert fasta.records[src].length == len(text)
        for _ in range(50):
            start = oracle.rng.randrange(len(text))
            end = oracle.rng.randrange(start, len(text) + 1)
            assert fasta.fetch(src, start, end) == text[start:end]
        assert fasta.fetch(src, 0, len(text)) == text
    fasta.close()


def test_no_newline_at_end(tmp_path):
    path = tmp_path / "x.fa"
    path.write_text(">a\nACGT\nAC")
    fasta = Fasta(str(path))
    assert fasta.fetch("a", 0, 6) == "ACGTAC" and fasta.fetch("a", 3, 5) == "TA"
    fasta.close()


def test_fai_is_used_only_if_current(tmp_path):
    path = tmp_path / "x.fa"
    path.write_text(">a\nACGT\n")
    index = tmp_path / "x.fa.fai"
    index.write_text("other\t4\t3\t4\t5\n")
    now = os.stat(path).st_mtime_ns
    for delta, name in ((1000, "other"), (-1000, "a")):
        os.utime(index, ns=(now + delta, now + delta))
        fasta = Fasta(str(path))
        assert list(fasta.records) == [name] and fasta.fetch(name, 1, 3) == "CG"
        fasta.close()


@pytest.mark.parametrize(
    "content, message",
    [(">a\nACGT\nAC\nACGT\n", "differ in length"), (">a\nAC\nACGT\n", "differ in length"),
     (">a\nAC\n\nAC\n", "differ in length"), (">a\nAC\n>a\nAC\n", "twice"),
     ("ACGT\n", "not a FASTA")],
)
def test_malformed_fasta(tmp_path, content, message):
    path = tmp_path / "x.fa"
    path.write_text(content)
    with pytest.raises(MafError, match=message):
        Fasta(str(path))


def test_names(tmp_path):
    a = write_fasta(tmp_path / "a.fa", {"chr1": "A", "hg.chr2": "C"})
    b = write_fasta(tmp_path / "b.fa", {"chr1": "G", "mm.chr1": "T", "NC_000001.1": "A"})
    genomes = Genomes()
    genomes.add(a, "hg")
    genomes.add(b)
    found = {s: genomes.find(s) for s in ("hg.chr1", "hg.chr2", "mm.chr1", "mm.chr2",
                                          "rn.chr1", "chr1", "NC_000001.1")}
    assert {s: f and f[0].fetch(f[1], 0, 1) for s, f in found.items()} == {
        "hg.chr1": "A", "hg.chr2": "C", "mm.chr1": "T", "mm.chr2": None, "rn.chr1": None,
        "chr1": "G", "NC_000001.1": "A",
    }
    genomes.close()


# --- checks -----------------------------------------------------------------

def test_true_alignments_have_no_problems(oracle, genomes):
    for _ in range(10):
        block = oracle.random_block(gaps=0.4)
        out, stats, unknown = problems(block, genomes)
        assert out == [] and stats["checked"] == 4 and not unknown


def test_output_of_window_and_merge_is_correct(oracle, genomes):
    for window in window_blocks([oracle.random_block(columns=200)], 40, 20):
        assert kinds(window, genomes) == []
    a, b = (Block([oracle.row(i, s + 5 * i * (s > 1000), "x" * 50) for i in range(3)])
            for s in (1000, 1050))
    (merged,) = merge_blocks([a, b], MergeOptions(max_distance=50))
    assert "N" in merged.rows[1].text and kinds(merged, genomes) == []


def test_lower_case_and_u(oracle, genomes):
    row = oracle.row(0, 100, "x" * 50)
    row.text = row.text.lower().replace("t", "u")
    assert kinds(Block([row]), genomes) == []


def test_soft_masked_genome(oracle, tmp_path):
    path = write_fasta(tmp_path / "sp0.fa", {"chr1": oracle.genomes["sp0.chr1"].lower()})
    genomes = Genomes()
    genomes.add(path, "sp0")
    for strand in "+-":
        assert kinds(Block([oracle.row(0, 100, "x" * 50, strand)]), genomes) == []
    genomes.close()


@pytest.mark.parametrize("strand", "+-")
def test_each_mistake_is_named(oracle, genomes, strand):
    def check(change):
        row = oracle.row(1, 700, "xx--" * 20, strand)
        change(row)
        ((problem),) = problems(Block([row]), genomes)[0]
        return problem.problem, problem.detail

    flipped = "-" if strand == "+" else "+"
    assert check(lambda r: setattr(r, "strand", flipped)) == (
        "strand", f"bases match strand {strand} at the same start")
    assert check(lambda r: setattr(r, "start", 5000 - r.end))[0] == "start"
    assert "other end" in check(lambda r: setattr(r, "start", 5000 - r.end))[1]
    assert check(lambda r: setattr(r, "text", reverse_complement(r.text)))[0] == "orientation"
    for shift in (-10, -1, 1, 10):
        assert check(lambda r: setattr(r, "start", r.start + shift)) == (
            "start", f"shifted by {shift:+d}, expected 700")
    assert check(lambda r: setattr(r, "start", r.start + 11))[0] == "bases"
    assert check(lambda r: setattr(r, "src_size", 4999)) == (
        "srcSize", "the sequence has 5000 bases, not 4999")
    assert check(lambda r: setattr(r, "start", 4990))[0] == "range"


def test_differences_are_counted(oracle, genomes):
    row = oracle.row(0, 100, "x" * 50)
    swap = str.maketrans("ACGT", "CATG")
    row.text = row.text[:10].translate(swap) + row.text[10:]
    (problem,) = problems(Block([row]), genomes)[0]
    assert (problem.problem, problem.detail) == ("bases", "10 of 50 differ")


def test_short_rows_are_not_explained(oracle, genomes):
    row = oracle.row(0, 100, "x" * 11)
    row.start += 1
    assert kinds(Block([row]), genomes) == ["bases"]
    row = oracle.row(0, 100, "x" * 12)
    row.start += 1
    assert kinds(Block([row]), genomes) == ["start"]


def test_unknown_sources_are_counted_not_reported(oracle, genomes):
    block = oracle.random_block()
    block.rows[1].src = "zz.chr1"
    out, stats, unknown = problems(block, genomes)
    assert out == [] and stats["checked"] == 3 and unknown == {"zz.chr1": 1}


def test_malformed_lines_do_not_end_the_check(oracle, genomes):
    good = oracle.random_block()
    lines = format_block(good).splitlines(keepends=True)
    lines[2] = lines[2].replace(" - ", " -1 ").replace(" + ", " 1 ")
    bad = oracle.row(0, 100, "x" * 50)
    bad.start += 3
    text = "".join(lines) + format_block(Block([bad])) + "a\ns x 0 9 + 9 AC\n"
    out, stats, _ = problems(text, genomes)
    assert [(p.line, p.problem) for p in out] == [(3, "format"), (7, "start"), (11, "format")]
    assert stats["checked"] == 3 + 1 + 1 + 1


# --- command line -----------------------------------------------------------

@pytest.fixture
def files(oracle, tmp_path):
    folder = tmp_path / "genomes"
    folder.mkdir()
    for src, text in oracle.genomes.items():
        write_fasta(folder / (src.split(".")[0] + ".fa"), {"chr1": text})
    maf = tmp_path / "x.maf"
    maf.write_text(format_block(oracle.random_block()))
    return maf, folder


def run(capsys, *argv):
    code = main([str(a) for a in argv])
    out, err = capsys.readouterr()
    return code, out, err


def test_cli(capsys, files):
    maf, folder = files
    code, out, err = run(capsys, "verify", "-d", folder, maf)
    assert code == 0 and out.count("\n") == 1 and "4 sequences checked, 0 with problems" in err
    code, out, err = run(capsys, "verify", "-g", f"sp1={folder / 'sp1.fa'}", maf)
    assert code == 0 and "1 sequences checked" in err and "no genome for 3 names" in err


def test_cli_problems_set_the_exit_status(capsys, files):
    maf, folder = files
    row = Row("sp1.chr1", 0, 4, "+", 4999, "ACGT")
    maf.write_text(format_block(Block([row])))
    code, out, err = run(capsys, "verify", "-d", folder, maf)
    assert code == 1 and "srcSize: 1" in err and out.splitlines()[1].startswith("1\tsp1.chr1\t0")


def test_cli_nothing_checked_is_an_error(capsys, files, tmp_path):
    maf, folder = files
    other = write_fasta(tmp_path / "o.fa", {"zz.chr1": "ACGT"})
    code, _, err = run(capsys, "verify", "-g", other, maf)
    assert code == 1 and "nothing was checked" in err
    with pytest.raises(SystemExit):
        main(["verify", str(maf)])


def test_cli_errors(capsys, files, tmp_path):
    maf, folder = files
    packed = tmp_path / "p.fa"
    packed.write_bytes(b"\x1f\x8b rest")
    code, _, err = run(capsys, "verify", "-g", packed, maf)
    assert code == 1 and "decompress" in err
    assert run(capsys, "verify", "-g", "sp1=/does/not/exist", maf)[0] == 1
    genome = folder / "sp1.fa"
    before = genome.read_text()
    assert run(capsys, "verify", "-d", folder, "-o", genome, maf)[0] == 1
    assert genome.read_text() == before
