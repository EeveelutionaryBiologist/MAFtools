# Changelog

## 2.0.0 (in development)

Rewritten as a package with a test suite. No dependencies.

### Wrong results in 1.0.0, now fixed

- `select`: all sequences received the start coordinate of the reference.
- `select`: the last column was lost when a region reached the end of a block.
- `select`: a block was lost if any sequence had only gaps in the region; regions of 1 bp were ignored.
- `select`: `chr1` also matched `chr10`, and did not match `hg38.chr1`.
- `select`: `--left-flank` and `--right-flank` had no effect.
- `select`, `merge`: a file header was written for every block.
- `filter`: `--min-length` kept only blocks of exactly that length.
- `filter`, `window`: the strand was written as `1` or `-1`; scores were replaced by 0.
- `window`: sequences without any base in the window were written with size 0.
- `merge`: could join sequences of different species.
- `merge`: discarded the first block if fewer than two sequences could be joined.
- `merge`: wrote wrong sizes if sequences differed in length.
- `merge`: joined overlapping and out-of-order blocks.
- `merge`: `--min-seqs` did not apply to merging, and shared sequences were counted twice.
- `merge`: `--max-length` ignored inserted bases.
- `merge`: crashed if the last block had too few sequences.
- `toBed`: coordinates of minus-strand references were not converted.
- `toBed`, `merge`: appended to existing output files.
- `describe`, `merge`: crashed on files without blocks.

### Added

- `species`: keep, drop or reorder species.
- `toFasta`: aligned or ungapped FASTA, per block or concatenated by species.
- `index`: random access for `select`.
- Gzipped input and output, and standard input.

### Changed

- Output goes to standard output unless `-o` is given, messages to standard error.
- BED names must match exactly, with or without the species prefix. Prefixes such as
  `NC_037638` for `NC_037638.1` no longer match.
- Malformed MAF and BED input is rejected with file name and line number.
- Windows without reference bases are skipped.
- Scores and `i`/`e` lines are kept for unchanged blocks and removed from cut or merged ones.
