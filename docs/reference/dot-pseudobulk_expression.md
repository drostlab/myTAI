# Create Pseudobulk Expression Data

Aggregate single-cell data by cell identity to create pseudobulk
expression.

## Usage

``` r
.pseudobulk_expression(expression, groups)
```

## Arguments

- expression:

  Expression matrix with genes as rows and cells as columns

- groups:

  Factor vector indicating which group each cell belongs to

## Value

Matrix of pseudobulked expression
