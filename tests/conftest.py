"""Oracle: synthetic genomes against which every alignment row can be checked."""

import random
from io import StringIO

import pytest

from maftools.blocks import Block, Row, count_bases
from maftools.mafio import MafReader, format_block

COMPLEMENT = str.maketrans("ACGT", "TGCA")


class Oracle:
    def __init__(self, seed=1, species=5, length=5000):
        self.rng = random.Random(seed)
        self.genomes = {
            f"sp{i}.chr1": "".join(self.rng.choices("ACGT", k=length)) for i in range(species)
        }

    def strand_sequence(self, src, strand):
        genome = self.genomes[src]
        return genome if strand == "+" else genome.translate(COMPLEMENT)[::-1]

    def row(self, i, start, mask, strand="+"):
        """Row of species i; mask is a string with '-' for gaps, anything else for bases."""
        src = f"sp{i}.chr1"
        bases = iter(self.strand_sequence(src, strand)[start:])
        text = "".join("-" if c == "-" else next(bases) for c in mask)
        return Row(src, start, count_bases(text), strand, len(self.genomes[src]), text)

    def random_block(self, start=1000, columns=100, rows=4, gaps=0.2, strands="+-"):
        out = []
        for i in range(rows):
            p = gaps / 2 if i == 0 else gaps
            mask = "".join("-" if self.rng.random() < p else "x" for _ in range(columns))
            strand = self.rng.choice(strands)
            out.append(self.row(i, start + 37 * i, mask, strand))
        return Block(out, {"score": "1234.0"})

    def check(self, block, bridged=False):
        """Assert that a block is valid MAF and true to the genomes."""
        assert block.rows
        (reread,) = MafReader(StringIO(format_block(block)))  # validates sizes and lengths
        assert reread.rows == block.rows
        for row in block.rows:
            assert row.size > 0
            bases = row.text.replace("-", "")
            truth = self.strand_sequence(row.src, row.strand)[row.start : row.end]
            assert len(truth) == len(bases)
            for got, want in zip(bases, truth):
                assert got == want or (bridged and got == "N"), (row.src, row.start)


@pytest.fixture(params=range(5))
def oracle(request):
    """Run with several seeds."""
    return Oracle(seed=request.param)
