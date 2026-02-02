# Downsample Expression Matrix by Groups

Downsample an expression matrix by randomly selecting a specified number
of samples from each group.

## Usage

``` r
downsample_expression(expression_matrix, groups, downsample = 10)
```

## Arguments

- expression_matrix:

  Expression matrix with genes as rows and samples as columns

- groups:

  Factor vector indicating which group each sample belongs to

- downsample:

  Integer, number of samples to keep per group (default: 10)

## Value

A dense expression matrix (genes x downsampled samples)

## Details

This function randomly samples up to `downsample` samples from each
group level. The returned expression matrix is converted to dense format
and maintains column names from the original matrix. Useful for creating
balanced subsets for visualization or when memory is limited.
