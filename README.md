# dekupl-computeNF

computeNF is a small utilitary that compute Normalization Factors (NF) from k-mer counts table. NF is the median of the ratios between sample counts and counts of a pseudo-reference obtained by taking the geometric mean of each k-mer across all samples.

# Installation

## Dependancies

- gcc compiler
- zlib

## Installation from the sources

1. Clone the project: `git clone https://github.com/Transipedia/dekupl-computeNF.git`
2. Compile the software: `make`
3. Place computeNF binary somewhere accessible to your `$PATH`

# Documentation

Usage:

    computeNF [-s SAMPLING_RATE] counts.tsv[.gz]

computeNF expect only one argument: the k-mers counts table (possibly gzipped). This table is a TSV (Tabulated-Separated Value) file, with an header line containing the sample names, and onr k-mer by row with its associated counts for each sample.

## Options

- `-s SAMPLING_RATE` option define a sampling rate at which k-mers are considered to compute NF. This will accelerate the processing time.

## Input example:

    tag     sample1 sample2 sample3 sample4
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA 627     357     345     446
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC 2       6       9       10
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG 16      15      25      28
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT 11      10      11      10
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAGA 10      9       20      14
    ...

## Output example

    sample  normalization_factor
    sample1 0.9831
    sample2 1.2342
    sample3 1.0023
