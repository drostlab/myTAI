# Count Transformation Functions

Predefined list of transformation functions for count data
normalization.

## Usage

``` r
COUNT_TRANSFORMS
```

## Format

A named list of transformation functions:

- none:

  Identity transformation (no change)

- sqrt:

  Square root transformation

- log2:

  log2(x+1) transformation

- vst:

  Variance stabilizing transformation (DESeq2)

- rlog:

  Regularized log transformation (DESeq2)

- rank:

  Rank transformation within each sample

## Details

This object provides a list of predefined transformation functions for
gene expression matrix. For \`rlog()\` and \`vst()\` transformations
(from \`DESeq2\`), the expression matrix is multiplied by 2 prior to
rounding. This is done to preserve more information for the lower
expression values (especially for TPM).
