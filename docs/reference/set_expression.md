# Gene Expression Transformation Functions

Collection of functions for transforming gene expression data in
PhyloExpressionSet objects.

Generic function to set the expression matrix in a PhyloExpressionSet
object.

## Usage

``` r
set_expression(phyex_set, new_expression, new_name = NULL, ...)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- new_expression:

  Matrix to set as the new expression

- new_name:

  Optional new name for the dataset

- ...:

  Additional arguments

## Value

A PhyloExpressionSet object with updated expression
