# Transform Expression Counts in PhyloExpressionSet

Apply a transformation function to the expression counts in a
PhyloExpressionSet.

## Usage

``` r
transform_counts(
  phyex_set,
  FUN,
  FUN_name = deparse(substitute(FUN)),
  new_name = NULL,
  ...
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- FUN:

  Function to apply

- FUN_name:

  Name of the transformation function (optional)

- new_name:

  Name for the new dataset (optional)

- ...:

  Additional arguments passed to FUN

## Value

A PhyloExpressionSet object with transformed expression data
