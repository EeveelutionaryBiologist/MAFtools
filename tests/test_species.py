from collections import Counter

import pytest

from maftools.blocks import Block
from maftools.operations import SpeciesOptions, species_blocks


def species(block, *names, stats=None, seen=None, **kw):
    return list(species_blocks([block], SpeciesOptions(names, **kw), stats, seen))


def sources(block):
    return [row.src for row in block.rows]


def test_keep_and_exclude_are_true_to_genomes(oracle):
    block = oracle.random_block(gaps=0.6)
    (kept,) = species(block, "sp0", "sp2.chr1")
    (left,) = species(block, "sp1", "sp3", exclude=True)
    assert sources(kept) == ["sp0.chr1", "sp2.chr1"] and kept == left
    oracle.check(kept)
    assert all(set(column) != {"-"} for column in zip(*(r.text for r in kept.rows)))


def test_names_match_whole_species_or_source(oracle):
    block = oracle.random_block()
    assert species(block, "sp") == species(block, "sp0.chr") == species(block, "chr1") == []


def test_unchanged_block_passes_through(oracle):
    block = oracle.random_block()
    (out,) = species(block, "sp3", "sp2", "sp1", "sp0")
    assert out is block


def test_reorder_makes_first_species_the_plus_strand_reference(oracle):
    for _ in range(10):
        block = oracle.random_block(gaps=0.3)
        (out,) = species(block, "sp2", "sp0", "sp3", reorder=True)
        assert sources(out) == ["sp2.chr1", "sp0.chr1", "sp3.chr1"]
        assert out.reference.strand == "+"
        oracle.check(out)
        assert out.reference.plus_interval() == block.rows[2].plus_interval()


def test_reorder_drops_blocks_without_the_new_reference(oracle):
    block = oracle.random_block(rows=3)
    assert species(block, "sp4", "sp0", reorder=True) == []


def test_require_all_and_min_seqs(oracle):
    block = oracle.random_block(rows=3)
    assert len(species(block, "sp0", "sp4")) == 1
    assert species(block, "sp0", "sp4", require_all=True) == []
    assert len(species(block, "sp0", "sp1", require_all=True)) == 1
    assert len(species(block, "sp0", "sp1", min_seqs=2)) == 1
    assert species(block, "sp0", "sp1", min_seqs=3) == []
    assert species(block, "sp0", "sp1", "sp2", exclude=True) == []


def test_stats_and_seen(oracle):
    stats, seen = Counter(), set()
    species(oracle.random_block(rows=3), "sp1", "sp9", stats=stats, seen=seen)
    assert stats == {"read": 1, "kept": 1, "changed_reference": 1} and seen == {"sp1"}


def test_reverse_complement(oracle):
    block = oracle.random_block(gaps=0.3)
    flipped = block.reverse_complement()
    oracle.check(flipped)
    assert [r.plus_interval() for r in flipped.rows] == [r.plus_interval() for r in block.rows]
    assert flipped.reverse_complement().rows == block.rows


def test_reverse_complement_keeps_case_and_ambiguity_codes(oracle):
    row = oracle.row(0, 0, "xxxx")
    row.text, row.quality = "aR-nC", "12-34"
    flipped = row.reverse_complement()
    assert (flipped.text, flipped.quality) == ("Gn-Yt", "43-21")
