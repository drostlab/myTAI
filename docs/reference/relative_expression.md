# Relative Expression Functions

Functions for computing and plotting relative expression profiles using
PhyloExpressionSet S7 objects.

Computes the relative expression profile for a given gene expression
matrix. The relative expression is calculated by normalizing the column
means of the matrix to a \[0, 1\] scale.

## Usage

``` r
relative_expression(count_matrix)
```

## Arguments

- count_matrix:

  A numeric matrix where columns represent developmental stages/cell
  types and rows represent genes.

## Value

A numeric vector of relative expression values for each stage (column)
in the input matrix.
