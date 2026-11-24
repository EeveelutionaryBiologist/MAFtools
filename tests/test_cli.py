import gzip
import subprocess
import sys
from pathlib import Path

import pytest

from maftools import __version__, cli
from maftools.cli import main
from maftools.mafio import MafReader

ROOT = Path(__file__).parent.parent
MAF, BED = str(ROOT / "Example/test.maf"), str(ROOT / "Example/test.bed")


def blocks(text):
    return list(MafReader(iter(text.splitlines(keepends=True))))


def run(capsys, *argv):
    code = main([str(a) for a in argv])
    out, err = capsys.readouterr()
    return code, out, err


def test_describe(capsys):
    assert run(capsys, "describe", MAF)[1] == (
        "Identifiers: 13\nSequences: 21\nNumber of blocks: 4\nColumns: 309\n"
        "Average sequences per block: 9.25\nAligned reference nucleotides: 297\n"
        "Average block length: 77.25\n"
    )


def test_blockdescribe(capsys):
    lines = run(capsys, "blockdescribe", MAF)[1].splitlines()
    assert lines[1] == "NC_037638.1\t1055913\t76\t78\t9\t0.6540\t0.1225" and len(lines) == 5


def test_to_bed_overwrites(capsys, tmp_path):
    out = tmp_path / "x.bed"
    for option in ("-o", "--output_bed"):
        run(capsys, "toBed", option, out, MAF)
    assert out.read_text().splitlines()[0] == "NC_037638.1\t1055913\t1055989"
    assert len(out.read_text().splitlines()) == 4


def test_to_bed_converts_minus_strand(capsys, tmp_path):
    maf = tmp_path / "m.maf"
    maf.write_text("a\ns x.1 10 4 - 100 ACGT\n")
    assert run(capsys, "toBed", maf)[1] == "x.1\t86\t90\n"


@pytest.mark.parametrize(
    "argv, count",
    [
        (["select", "-b", BED], 4),
        (["filter", "-g", "1", "-m", "0", "-l", "79", "-n", "1"], 2),
        (["window", "-l", "40", "-s", "20"], None),
        (["merge", "-d", "100"], None),
    ],
)
def test_maf_output(capsys, argv, count):
    code, out, _ = run(capsys, *argv, MAF)
    assert code == 0 and out.count("##maf") == 1
    assert len(blocks(out)) == (count or len(blocks(out))) > 0


def test_commands_can_be_chained(capsys, tmp_path):
    step = tmp_path / "w.maf"
    run(capsys, "window", "-l", "40", "-s", "40", "-o", step, MAF)
    code, out, err = run(capsys, "merge", step)
    assert code == 0 and blocks(out)


def test_input_is_never_overwritten(capsys, tmp_path):
    maf = tmp_path / "m.maf"
    maf.write_text(Path(MAF).read_text())
    code, _, err = run(capsys, "filter", "-o", maf, maf)
    assert code == 1 and "input" in err and maf.read_text() == Path(MAF).read_text()


@pytest.mark.parametrize("command", ["describe", "blockdescribe", "toBed", "merge", "window"])
def test_empty_file(capsys, tmp_path, command):
    maf = tmp_path / "e.maf"
    maf.write_text("##maf version=1\n\n")
    assert run(capsys, command, maf)[0] == 0


def test_errors_are_messages(capsys, tmp_path):
    binary = tmp_path / "b.maf"
    binary.write_bytes(b"\x1f\x8b\x08\xff\xfe")
    for path in (tmp_path / "missing.maf", binary):
        code, _, err = run(capsys, "describe", path)
        assert code == 1 and err.startswith("error:")


@pytest.mark.parametrize(
    "argv",
    [[], ["nonsense", MAF], ["select", MAF], ["filter", "-g", "2", MAF],
     ["window", "-l", "10", "-s", "20", MAF], ["window", "-s", "0", MAF]],
)
def test_usage_errors(capsys, argv):
    with pytest.raises(SystemExit) as error:
        main(argv)
    assert error.value.code == 2


def test_version(capsys):
    with pytest.raises(SystemExit):
        main(["--version"])
    assert __version__ in capsys.readouterr().out


def test_script_runs_from_checkout():
    done = subprocess.run([sys.executable, ROOT / "MAFtools.py", "describe", MAF],
                          capture_output=True, text=True)
    assert done.returncode == 0 and "Number of blocks: 4" in done.stdout


# --- compression, pipes, index ----------------------------------------------

def shell(*argv, stdin=b""):
    return subprocess.run([sys.executable, ROOT / "MAFtools.py", *map(str, argv)],
                          input=stdin, capture_output=True)


def test_gzip_input_and_output(capsys, tmp_path):
    packed = tmp_path / "any_name"
    packed.write_bytes(gzip.compress(Path(MAF).read_bytes()))
    plain = run(capsys, "describe", MAF)[1]
    assert run(capsys, "describe", packed)[1] == plain
    out = tmp_path / "o.maf.gz"
    run(capsys, "window", "-o", out, packed)
    assert gzip.decompress(out.read_bytes()).decode() == run(capsys, "window", MAF)[1]
    assert out.read_bytes()[4:8] == bytes(4)  # no time stamp: output is reproducible


def test_multi_member_gzip(capsys, tmp_path):
    data = Path(MAF).read_bytes()
    packed = tmp_path / "x.maf.gz"
    packed.write_bytes(gzip.compress(data[:3000]) + gzip.compress(data[3000:]))
    assert run(capsys, "describe", packed)[1] == run(capsys, "describe", MAF)[1]


def test_damaged_gzip_is_an_error(capsys, tmp_path):
    packed = tmp_path / "x.maf.gz"
    data = gzip.compress(Path(MAF).read_bytes())
    for bad in (data[: len(data) // 2], data[:200] + bytes(100) + data[300:]):
        packed.write_bytes(bad)
        code, _, err = run(capsys, "describe", packed)
        assert code == 1 and "error" in err and "Traceback" not in err


@pytest.mark.parametrize("pack", [bytes, gzip.compress])
def test_standard_input(capsys, pack):
    done = shell("window", "-l", "40", "-", stdin=pack(Path(MAF).read_bytes()))
    assert done.returncode == 0
    assert done.stdout.decode() == run(capsys, "window", "-l", "40", MAF)[1]


def test_bed_from_standard_input(capsys):
    done = shell("select", "-b", "-", MAF, stdin=Path(BED).read_bytes())
    assert done.stdout.decode() == run(capsys, "select", "-b", BED, MAF)[1]
    assert shell("select", "-b", "-", "-").returncode == 2


def test_closed_pipe_is_quiet():
    command = f"'{sys.executable}' '{ROOT / 'MAFtools.py'}' window -l 5 -s 1 '{MAF}' | head -1"
    done = subprocess.run(command, shell=True, capture_output=True, text=True)
    assert done.stdout == "##maf version=1\n" and done.stderr == ""


def test_index_is_used_and_checked(capsys, tmp_path, monkeypatch):
    maf = tmp_path / "x.maf"
    maf.write_text(Path(MAF).read_text())
    plain = run(capsys, "select", "-l", "5", "-b", BED, maf)
    assert run(capsys, "index", maf)[0] == 0
    assert run(capsys, "select", "-l", "5", "-b", BED, maf) == plain

    calls = []
    original = cli.MafIndex.blocks
    monkeypatch.setattr(cli.MafIndex, "blocks",
                        lambda *a: calls.append(1) or original(*a))
    run(capsys, "select", "-b", BED, maf)
    assert calls == [1]

    maf.write_text(maf.read_text() + "\n")
    code, out, err = run(capsys, "select", "-l", "5", "-b", BED, maf)
    assert out == plain[1] and "out of date" in err and calls == [1]


def test_flanks_reach_blocks_through_the_index(capsys, tmp_path):
    maf, bed = tmp_path / "x.maf", tmp_path / "x.bed"
    maf.write_text(Path(MAF).read_text())
    bed.write_text("NC_037638.1\t2006350\t2006360\n")  # 3 bases after a block
    run(capsys, "index", maf)
    assert len(blocks(run(capsys, "select", "-l", "5", "-b", bed, maf)[1])) == 1


def test_damaged_index_is_ignored(capsys, tmp_path):
    maf = tmp_path / "x.maf"
    maf.write_text(Path(MAF).read_text())
    Path(str(maf) + ".mfi").write_text("not a database, but long enough to look like one" * 50)
    code, out, err = run(capsys, "select", "-b", BED, maf)
    assert code == 0 and len(blocks(out)) == 4 and "not used" in err


# --- species and FASTA ------------------------------------------------------

def test_species(capsys, tmp_path):
    names = tmp_path / "names.txt"
    names.write_text("# keep\nBombus_terrestris\n\n")
    code, out, err = run(capsys, "species", "-s", "NC_037638", "-f", names, "-s", "Nobody", MAF)
    assert code == 0 and "not found: Nobody" in err
    assert {row.species for b in blocks(out) for row in b.rows} == {"NC_037638",
                                                                    "Bombus_terrestris"}
    assert "lost their reference" in run(capsys, "species", "-s", "Bombus_terrestris", MAF)[2]
    reordered = blocks(run(capsys, "species", "-r", "-s", "Bombus_terrestris,NC_037638", MAF)[1])
    assert {(b.reference.species, b.reference.strand) for b in reordered} == {
        ("Bombus_terrestris", "+")}
    assert run(capsys, "species", MAF)[0] == 1
    with pytest.raises(SystemExit):
        main(["species", "-x", "-r", "-s", "a", MAF])


def test_list_is_not_overwritten(capsys, tmp_path):
    names = tmp_path / "names.txt"
    names.write_text("NC_037638\n")
    assert run(capsys, "species", "-f", names, "-o", names, MAF)[0] == 1
    assert names.read_text() == "NC_037638\n"


def test_to_fasta(capsys):
    assert run(capsys, "toFasta", MAF)[1].count(">") == 37
    out = run(capsys, "toFasta", "-c", MAF)[1]
    assert out.count(">") == 13
    ungapped = run(capsys, "toFasta", "-u", MAF)[1].splitlines()
    assert not any("-" in line for line in ungapped if not line.startswith(">"))
    lengths = {sum(map(len, r.split("\n")[1:])) for r in out.split(">")[1:]}
    assert lengths == {309}


def test_dedup(capsys, tmp_path):
    maf = tmp_path / "d.maf"
    maf.write_text("a\ns x.1 0 4 + 9 ACGT\ns y.1 0 4 + 9 TTTT\ns y.2 0 4 + 9 ACGA\n")
    code, out, err = run(capsys, "dedup", maf)
    assert code == 0 and [r.src for r in blocks(out)[0].rows] == ["x.1", "y.2"]
    assert err == ("1 of 1 blocks had duplicates; 1 sequences removed, 0 blocks dropped\n")
    assert [r.src for r in blocks(run(capsys, "dedup", "-k", "first", maf)[1])[0].rows] == [
        "x.1", "y.1"]
    assert len(blocks(run(capsys, "dedup", "-k", "none", maf)[1])[0].rows) == 1
    assert blocks(run(capsys, "dedup", "-n", "3", maf)[1]) == []
    assert "see dedup" in run(capsys, "toFasta", "-c", maf)[2]
