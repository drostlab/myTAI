# Plot Mean Relative Expression Levels as Barplot

Plots mean relative expression levels for age category groups using a
PhyloExpressionSet S7 object, with statistical testing.

## Usage

``` r
plot_relative_expression_bar(phyex_set, groups, p_adjust_method = NULL, ...)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet).

- groups:

  A list of integer vectors specifying age categories (e.g.,
  phylostrata) for each group (2+ groups).

- p_adjust_method:

  P-value adjustment for multiple testing.

- ...:

  Further arguments passed to ggplot2 geoms.

## Value

ggplot2 object.
