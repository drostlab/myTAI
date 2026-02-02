# QQ plot comparing partial TAI distributions across developmental stages against a reference stage

*plot_distribution_partialTAI_qqplot* generates a QQ plot to compare the
partial TAI distributions of various developmental stages against a
reference stage. It visualizes quantile differences between the
reference and other stages, highlights each stage with distinct colors,
and annotates the plot with the p-values from the nonparametric
[`ks.test`](https://rdrr.io/r/stats/ks.test.html) to indicate the
significance of distribution differences.

## Usage

``` r
plot_distribution_pTAI_qqplot(
  phyex_set,
  reference_stage_index = 1,
  xlab = "Quantiles of Reference Stage",
  ylab = "Quantiles of Other Stages",
  main = "QQ Plot: Developmental Stages vs Reference Stage (p-values from KS test)",
  alpha = 0.7,
  size = 1.2
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet).

- reference_stage_index:

  An integer specifying the index of the reference developmental stage.
  The partial TAI distribution of this stage will be used as the
  reference for comparisons with other stages (default: stage index 1).

- xlab:

  Label of x-axis.

- ylab:

  Label of y-axis.

- main:

  Figure title.

- alpha:

  Transparency of the points.

- size:

  Size of the points.

## Value

A ggplot2 object showing a qqplot of partial TAI distributions

## Author

Filipa Martins Costa
