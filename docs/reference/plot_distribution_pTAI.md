# Partial TAI Distribution Plotting Functions

Functions for plotting and comparing partial TAI distributions using
PhyloExpressionSet S7 objects.

*plot_distribution_pTAI* generates 2 plots that help to compare the
distribution of the quotient of expression by partial TAI through
various developmental stages or cell types, highlighting each stage with
distinct colors.

## Usage

``` r
plot_distribution_pTAI(
  phyex_set,
  stages = NULL,
  xlab = "Expression / Partial TAI",
  ylab = "Density",
  main = "Density Distribution of Expression / Partial TAI by Developmental Stage"
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet).

- stages:

  A numeric vector specifying the indices of the stages to compare. Each
  index corresponds to a stage in the PhyloExpressionSet. If NULL, all
  stages are used.

- xlab:

  Label of x-axis.

- ylab:

  Label of y-axis.

- main:

  Figure title.

## Value

A ggplot2 object showing partial TAI distributions

## Author

Filipa Martins Costa
