from io import StringIO

import pytest

from maftools.mafio import MafError, MafReader, format_block, write_maf

FULL = """##maf version=1 scoring=multiz
# a tree

a score=5.0 pass=2
s hg.chr1 10 4 + 100 AC-GT
q hg.chr1 99-99
i hg.chr1 C 0 C 0
s mm.chr2 20 5 - 200 ACCGT
e rn.chr3 5 10 + 300 I

"""


def read(text):
    return list(MafReader(StringIO(text), "x.maf"))


def test_roundtrip_keeps_everything():
    reader = MafReader(StringIO(FULL))
    out = StringIO()
    assert write_maf(out, reader, reader.header) == 1
    assert [l.split() for l in out.getvalue().split("\n")] == [l.split() for l in FULL.split("\n")]


def test_tabs_missing_header_and_missing_blank_line():
    blocks = read("\na\ns x.1 \t0 \t2 \t+ \t9 \tAC\na\ns x.1 2 2 + 9 GT")
    assert [b.reference.start for b in blocks] == [0, 2]


def test_empty_file():
    assert read("") == [] and read("##maf version=1\n\n") == []


def test_header_is_added_if_missing():
    out = StringIO()
    write_maf(out, [], ["# only a comment"])
    assert out.getvalue().startswith("##maf version=1\n# only a comment\n")


@pytest.mark.parametrize(
    "body, message",
    [
        ("a\ns x 0 2 1 9 AC\n", "strand"),
        ("a\ns x 0 3 + 9 AC\n", "3 but the sequence has 2"),
        ("a\ns x 0 2 + 9 AC\ns y 0 3 + 9 ACG\n", "differ in length"),
        ("a\ns x 0 2 + 9\n", "expected 7"),
        ("a\ns x zero 2 + 9 AC\n", "integers"),
        ("a\ns x 0 2 + 9 AC\ni y C 0 C 0\n", "does not follow"),
        (">fasta\nACGT\n", "expected an alignment block"),
    ],
)
def test_malformed_input_is_reported_with_position(body, message):
    with pytest.raises(MafError, match=message) as error:
        read(body)
    assert str(error.value).startswith("x.maf:")


def test_output_is_readable_by_biopython(oracle):
    AlignIO = pytest.importorskip("Bio.AlignIO")
    block = oracle.random_block()
    text = "##maf version=1\n\n" + format_block(block)
    (alignment,) = AlignIO.parse(StringIO(text), "maf")
    for record, row in zip(alignment, block.rows, strict=True):
        a = record.annotations
        assert (record.id, a["start"], a["size"]) == (row.src, row.start, row.size)
        assert a["strand"] == (1 if row.strand == "+" else -1)
