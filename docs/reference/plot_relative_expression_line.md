# Plot Relative Expression Profiles (Line Plot)

Plots relative expression profiles for age categories using a
PhyloExpressionSet S7 object.

## Usage

``` r
plot_relative_expression_line(
  phyex_set,
  groups,
  modules = NULL,
  adjust_range = TRUE,
  alpha = 0.1,
  ...
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet).

- groups:

  A list of integer vectors specifying age categories (e.g.,
  phylostrata) for each group (1 or 2 groups).

- modules:

  Optional list for shading modules: list(early=..., mid=..., late=...).

- adjust_range:

  Logical, adjust y-axis range for both panels (if 2 groups).

- alpha:

  Transparency for shaded module area.

- ...:

  Further arguments passed to ggplot2 geoms.

## Value

ggplot2 object or list of ggplot2 objects.
