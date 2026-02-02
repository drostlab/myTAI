# Comparing expression levels distributions across developmental stages

`plot_distribution_expression` generates plots that help to compare the
distribution of expression levels through various developmental stages
or cell types, highlighting each stage with distinct colors. By default,
a log transformation is applied to the expression values.

## Usage

``` r
plot_distribution_expression(
  phyex_set,
  show_identities = TRUE,
  show_strata = FALSE,
  log_transform = TRUE,
  seed = 123
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet).

- show_identities:

  Logical, whether to show identity-specific distributions.

- show_strata:

  Logical, whether to show stratum-specific distributions.

- log_transform:

  Logical, whether to apply log transformation to expression values
  (default: TRUE).

- seed:

  Seed for reproducible color selection.

## Value

A ggplot2 object showing expression levels distributions across
identities

## Recommendation

Apply a square root transformation to enhance the visualization of
differences in the distributions:
`plot_distribution_expression(transform_counts(phyex_set, sqrt))`

## Author

Filipa Martins Costa
