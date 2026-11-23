# MAFtools

A comprehensive Python toolkit for processing and analyzing Multiple Alignment Format (MAF) files.

## Overview

MAFtools is a command-line utility that provides various operations for MAF alignment files, including filtering, merging, windowing, and format conversion. It's particularly useful for comparative genomics analyses where you need to manipulate whole-genome alignments.

**Version:** 2.0.0 (in development)  
**Author:** Christopher Klapproth  
**Institution:** University Leipzig  
**License:** GPLv2

### Installation
Requires Python 3.10 or newer. There are no other dependencies.

```bash
pip install .
```

Without installing, use `python MAFtools.py` from a checkout in place of `maftools`.

> **Results from version 1.0.0:** `select`, `filter`, `window` and `merge` had bugs that
> produced wrong output, mostly without an error. See [CHANGELOG.md](CHANGELOG.md)
> before relying on files made with that version.

## Usage

```bash
maftools [program] [options] [MAF file]
```

Where `program` is one of: `describe`, `merge`, `filter`, `window`, `select`, `toBed`, `blockdescribe`

Results go to standard output unless `-o` is given; an existing output file is replaced.
Progress messages go to standard error. The first sequence of each block is the reference.

## Core Functionality

### 1. describe
Provides summary statistics about a MAF file.

```bash
maftools describe alignment.maf
```

**Output includes:**
- Number of unique identifiers and sequences
- Total number of alignment blocks
- Total columns and average block length
- Average sequences per block
- Aligned reference nucleotides

### 2. select
Extracts alignment blocks based on BED coordinates.

```bash
maftools select -b regions.bed -o extracted.maf alignment.maf
```

**Options:**
- `-b, --bed`: BED file with coordinates to extract (required)
- `-r, --right-flank`: Add nucleotides to the right side of each coordinate (default: 0)
- `-l, --left-flank`: Add nucleotides to the left side of each coordinate (default: 0)
- `-o, --output`: Output MAF file (default: standard output)

**Example:** If a BED entry is `species1 1560 1960` and a MAF block starts at position 1500 with length 500, the program extracts the overlapping region (positions 1560-1960).

**Behavior:**
- BED names match the full sequence name (`hg38.chr1`) or the part after the species (`chr1`).
  Versioned accessions such as `NC_037638.1` must be given in full.
- Sequences without any base in the region are left out.
- Insertions relative to the reference stay with the preceding reference base.
- The BED strand column is ignored; flanks refer to the plus strand.

### 3. filter
Filters MAF alignment blocks based on various criteria.

```bash
maftools filter -g 0.25 -n 3 -o filtered.maf alignment.maf
```

**Options:**
- `-g, --max-gaps`: Maximum gap fraction allowed (default: 0.25)
- `-a, --max-similarity`: Maximum similarity of any sequence to the reference (default: 1.0)
- `-m, --min-similarity`: Minimum similarity of every sequence to the reference (default: 0.5)
- `-l, --min-length`: Minimum block length in columns (default: 40)
- `-n, --min-seqs`: Minimum number of sequences required (default: 3)
- `-o, --output`: Output MAF file (default: standard output)

**Filter criteria:**
- Average gap fraction across all sequences
- Pairwise sequence similarity to the reference
- Number of sequences in the block
- Alignment block length

### 4. window
Applies a sliding window approach to create overlapping smaller alignment blocks.

```bash
maftools window -l 120 -s 40 -o windowed.maf alignment.maf
```

**Options:**
- `-l, --length`: Target length of cut alignment blocks (default: 120)
- `-s, --step`: Step size for sliding window (default: 40)
- `-o, --output`: Output MAF file (default: standard output)

**Behavior:**
- Blocks shorter than the window length are kept as-is
- Sequences without any base in a window are left out; windows without reference bases are skipped
- Longer blocks are cut into overlapping windows
- Example: A 200bp block with step=40 creates windows starting at positions 0, 40, 80, etc. The last window is shifted to end with the block.

### 5. merge
Merges adjacent alignment blocks based on species consensus and distance criteria.

```bash
maftools merge -s 0.75 -d 100 -o merged.maf alignment.maf
```

**Options:**
- `-s, --species-consensus`: Minimum consensus between blocks for merging (default: 0.75)
- `-r, --no-reference`: Apply `--max-gaps` to the first sequence as well
- `-d, --max-distance`: Maximum number of unaligned bases between blocks (default: 0)
- `-l, --max-length`: Maximum number of columns of merged blocks (default: 1000)
- `-g, --max-gaps`: Maximum gap fraction before dropping sequences (default: 0.5)
- `-m, --min-seqs`: Minimum number of sequences to merge, and to write a block (default: 2)
- `-o, --output`: Output MAF file (default: standard output)

**Behavior:**
- Sequences are joined only if name and strand agree and the distance is within `--max-distance`;
  others are dropped from the merged block.
- Unaligned bases between blocks are written as `N`, padded with gaps.

### 6. toBed
Converts MAF alignment blocks to BED format.

```bash
maftools toBed -o alignment.bed alignment.maf
```

**Options:**
- `-o, --output`: Output BED file (default: standard output)

**Output format:**
Each line represents one alignment block's coordinates:
```
sequence_name    start    end
```
Coordinates are always on the plus strand.

### 7. blockdescribe
Generates a detailed tab-delimited summary of each alignment block.

```bash
maftools blockdescribe -o blocks_summary.tsv alignment.maf
```

**Options:**
- `-o, --output`: Output TSV file (default: standard output)

**Output columns:**
- `refseq`: Reference sequence name (first sequence in block)
- `start`: Start coordinate of reference sequence
- `size`: Number of non-gap bases in reference sequence
- `length`: Total number of columns in the alignment block
- `seqs`: Number of sequences in the block
- `mean_pairwise_identity`: Average pairwise identity of all sequences to reference (0-1)
- `gaps_fraction`: Average fraction of gaps across all sequences (0-1)

**Example output:**
```tsv
refseq          start   size    length  seqs  mean_pairwise_identity  gaps_fraction
hg38.chr1       1000    450     500     5     0.8523                  0.1000
hg38.chr1       2000    380     400     7     0.9125                  0.0500
hg38.chr2       5000    290     300     4     0.7856                  0.0333
```

This is useful for:
- Quality assessment of alignment blocks
- Identifying highly conserved or divergent regions
- Pre-filtering analysis to select blocks of interest
- Generating summary statistics for downstream analysis



## Examples

### Extract specific genomic regions:

```
# Extract these regions from MAF
maftools select -b Example/test.bed -o extracted.maf Example/test.maf
```

### Filter alignments for quality or similarity criteria:
```bash
# Keep blocks with <15% gaps, 3+ sequences, and 70-95% similarity
maftools filter -g 0.15 -n 3 -m 0.70 -a 0.95 -o filtered.maf Example/test.maf
```

### Create sliding windows for analysis:
```bash
# Create 100bp windows, each starting 25bp after the previous one
maftools window -l 100 -s 25 -o windows.maf genome.maf
```

## MAF Format Notes

- Coordinates are 0-based, half-open intervals
- The reference sequence is typically the first sequence in each block
- Gap characters are represented by '-'
- For sequences on the minus strand, start counts from the end of the plus strand
- Scores and `i`/`e` lines are kept for unchanged blocks and removed from cut or merged ones
- Each block contains aligned sequences with metadata (start, size, strand, srcSize)

## Python API

```python
from maftools import MafReader

with open("alignment.maf") as handle:
    for block in MafReader(handle):
        print(block.reference.src, block.reference.start, block.columns)
```

`maftools.operations` has the same transformations as the command line, as generators.

## Development

```bash
pip install -e '.[test]'
pytest
```

## Citation

If you use MAFtools in your research, please cite:
```
Klapproth, C. (2024). MAFtools: A toolkit for Multiple Alignment Format file manipulation. 
University Leipzig. Version 1.0.0.
```
