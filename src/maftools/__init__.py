"""Filter, cut, merge and summarise Multiple Alignment Format (MAF) files."""

__version__ = "2.0.0.dev0"
__author__ = "Christopher Klapproth"

from .blocks import Block, Row
from .mafio import MafError, MafReader, format_block, write_maf

__all__ = ["Block", "MafError", "MafReader", "Row", "format_block", "write_maf"]
