# dekupl-computeNF

computeNF is a small utilitary that compute Normalization Factors (NF) from k-mer counts table. NF is the median of the ratios between sample counts and counts of a pseudo-reference obtained by taking the geometric mean of each k-mer across all samples.

In this implementation we compute the median using an incremental median estimator as proposed [here]( https://stackoverflow.com/a/2144754/3730728).

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

    computeNF counts.tsv[.gz]

computeNF expect only one argument: the k-mers counts table (possibly gzipped). This table is a TSV (Tabulated-Separated Value) file, with an header line containing the sample names, and onr k-mer by row with its associated counts for each sample.

## Input example:

    kmer  sample1 sample2 sample3
    AAG 12  23  45
    CGA 0 3 5
    ...

## Output example

    sample  normalization_factor
    sample1 0.9831
    sample2 1.2342
    sample3 1.0023


