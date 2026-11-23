import subprocess
import sys
from pathlib import Path

import pytest

from maftools import __version__
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
