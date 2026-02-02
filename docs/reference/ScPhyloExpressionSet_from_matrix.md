# Create Single-Cell PhyloExpressionSet from Expression Matrix

Create a ScPhyloExpressionSet object from an expression matrix and
metadata.

## Usage

``` r
ScPhyloExpressionSet_from_matrix(
  expression_matrix,
  strata,
  metadata,
  groups_column = NULL,
  name = "Single-Cell Phylo Expression Set",
  ...
)
```

## Arguments

- expression_matrix:

  Sparse or dense expression matrix with genes as rows and cells as
  columns

- strata:

  Factor vector of phylostratum assignments for each gene

- metadata:

  Data frame with cell metadata, rownames should match colnames of
  expression_matrix

- groups_column:

  Character string specifying which metadata column to use for initial
  grouping (default: first factor column found)

- name:

  A character string naming the dataset (default: "Single-Cell Phylo
  Expression Set")

- ...:

  Additional arguments passed to ScPhyloExpressionSet constructor

## Value

A ScPhyloExpressionSet object

## Details

This function creates a ScPhyloExpressionSet from basic components. The
`groups_column` parameter determines the initial `selected_idents`
value, which can be changed later using the setter. All discrete columns
in metadata are automatically converted to factors for consistent
handling.
