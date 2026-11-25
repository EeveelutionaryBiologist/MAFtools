from collections import Counter
from itertools import combinations

import pytest

from maftools.blocks import Block, Row
from maftools.cli import main
from maftools.mafio import format_block
from maftools.metrics import pairwise_identity
from maftools.operations import ThinOptions, thin_blocks, thin_rows


def rows_with_identities(*identities, length=100):
    """A reference of A, and rows that differ from it in their own columns."""
    rows = [Row("s0.1", 0, length, "+", 1000, "A" * length)]
    for i, identity in enumerate(identities, 1):
        changed = round(length * (1 - identity))
        rows.append(Row(f"s{i}.1", 0, length, "+", 1000,
                        "ACGT"[i % 3 + 1] * changed + "A" * (length - changed)))
    return rows


def thin(rows, **kw):
    return thin_rows(rows, ThinOptions(**{"max_identity": 1.0, **kw}))


def test_two_sequences_are_the_reference_and_the_closest_to_the_target():
    rows = rows_with_identities(0.5, 0.95, 0.82, 0.7)
    assert thin(rows, max_seqs=2) == [0, 3]
    assert thin(rows, max_seqs=2, target_identity=0.9) == [0, 2]
    assert thin(rows, max_seqs=2, target_identity=0.5) == [0, 1]


def test_every_step_is_the_best_possible(oracle):
    """Compare with a direct implementation of the definition."""
    rng = oracle.rng
    for _ in range(20):
        texts = ["".join(rng.choices("ACGT-", (6, 1, 1, 1, 2), k=40)) for _ in range(9)]
        rows = [Row(f"s{i}.1", 0, 0, "+", 99, text) for i, text in enumerate(texts)]
        selected = [0]
        while len(selected) < 5:
            def cost(i):
                pairs = list(combinations([*selected, i], 2))
                return sum((pairwise_identity(texts[a], texts[b]) - 0.8) ** 2
                           for a, b in pairs) / len(pairs)
            left = [i for i in range(9) if i not in selected]
            best = min(round(cost(i), 9) for i in left)
            selected.append(next(i for i in left if round(cost(i), 9) == best))
        assert thin(rows, max_seqs=5) == sorted(selected)


def test_near_identical_rows_are_dropped_even_below_the_limit():
    rows = rows_with_identities(1.0, 0.995, 0.99, 0.6, length=1000)
    assert thin_rows(rows, ThinOptions()) == [0, 3, 4]
    assert thin_rows(rows, ThinOptions(max_identity=1.0)) == [0, 1, 2, 3, 4]


def test_near_identical_rows_are_compared_with_each_other():
    rows = rows_with_identities(0.6, 0.6, 0.6)
    rows[2].text = rows[1].text
    assert thin_rows(rows, ThinOptions()) == [0, 1, 3]


def test_dropped_rows_do_not_drop_others():
    # b is close to a, c is close to b but not to a
    rows = rows_with_identities(0.95, 0.90)
    rows[1].text = "C" * 5 + "A" * 95
    rows[2].text = "C" * 10 + "A" * 90
    assert thin(rows, max_identity=0.93) == [0, 2]


def test_limits():
    rows = rows_with_identities(0.5, 0.95, 0.82, 0.7)
    assert thin(rows, max_seqs=1) == [0]
    assert thin(rows, max_seqs=5) == thin(rows, max_seqs=50) == [0, 1, 2, 3, 4]


def test_ties_go_to_the_first_row():
    rows = rows_with_identities(0.8, 0.8)
    assert thin(rows, max_seqs=2) == [0, 1]


def test_blocks_are_true_to_genomes(oracle):
    block = oracle.random_block(rows=5, gaps=0.5)
    stats = Counter()
    (out,) = thin_blocks([block], ThinOptions(max_seqs=3), stats)
    oracle.check(out)
    assert len(out.rows) == 3 and out.reference.src == block.reference.src
    assert [r.src for r in out.rows] == sorted(r.src for r in out.rows)
    assert all(set(column) != {"-"} for column in zip(*(r.text for r in out.rows)))
    assert stats == {"read": 1, "thinned": 1, "removed": 2}


def test_small_blocks_are_unchanged(oracle):
    block = oracle.random_block()
    (out,) = thin_blocks([block], ThinOptions())
    assert out is block


@pytest.mark.parametrize("pair, identity", [
    (("AC-G", "ac-T"), 2 / 3), (("--", "--"), 0), (("A-", "-A"), 0), (("", ""), 0),
    (("AC--", "ACGT"), 0.5),
])
def test_identity(pair, identity):
    assert pairwise_identity(*pair) == pytest.approx(identity)


def test_cli(capsys, tmp_path):
    maf = tmp_path / "x.maf"
    rows = rows_with_identities(0.5, 0.95, 0.82, 1.0)
    maf.write_text(format_block(Block(rows)))
    assert main(["thin", "-n", "2", str(maf)]) == 0
    out, err = capsys.readouterr()
    assert out.count("\ns ") == 2 and "s3.1" in out
    assert err == "1 of 1 blocks thinned; 3 sequences removed\n"
    main(["thin", "-n", "2", "-t", "0.5", str(maf)])
    assert "s1.1" in capsys.readouterr().out
    main(["thin", "-x", "1", str(maf)])
    assert capsys.readouterr().out.count("\ns ") == 5
