# Convert Data to BulkPhyloExpressionSet

Convert a data frame with phylostratum, gene ID, and expression data
into a BulkPhyloExpressionSet object.

## Usage

``` r
BulkPhyloExpressionSet_from_df(
  data,
  groups = colnames(data[, 3:ncol(data)]),
  name = "PhyloExpressionSet",
  strata_legend = NULL,
  ...
)
```

## Arguments

- data:

  A data frame where column 1 contains phylostratum information, column
  2 contains gene IDs, and columns 3+ contain expression data

- groups:

  A factor or character vector indicating which group each sample
  belongs to. Default uses column names from expression data

- name:

  A character string naming the dataset. Default uses the variable name

- strata_legend:

  A data frame with two columns: phylostratum assignments and name of
  each stratum. If NULL, no labels will be added (default: NULL) If
  NULL, uses sorted unique values from column 1

- ...:

  Additional arguments passed to BulkPhyloExpressionSet constructor

## Value

A BulkPhyloExpressionSet object
