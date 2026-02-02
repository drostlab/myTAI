# as_BulkPhyloExpressionSet

This function is an alias for `BulkPhyloExpressionSet_from_df`. Please
refer to the documentation for `BulkPhyloExpressionSet_from_df` for
usage details, arguments, and examples.

## Usage

``` r
as_BulkPhyloExpressionSet(
  data,
  groups = colnames(data[, 3:ncol(data)]),
  name = "PhyloExpressionSet",
  strata_legend = NULL,
  ...
)
```

## Arguments

- data:

  A data frame with phylostratum assignments and gene expression data

- groups:

  Vector of group labels for the samples/replicates

- name:

  Character string to name the dataset

- strata_legend:

  Optional data frame mapping phylostratum numbers to labels

- ...:

  Additional arguments passed to BulkPhyloExpressionSet_from_df

## Value

A BulkPhyloExpressionSet object

## See also

[`BulkPhyloExpressionSet_from_df`](https://drostlab.github.io/myTAI/reference/BulkPhyloExpressionSet_from_df.md)
