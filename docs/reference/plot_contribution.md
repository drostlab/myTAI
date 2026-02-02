# Plot Phylostratum Contribution to Transcriptomic Index

Create a stacked area plot showing the contribution of each phylostratum
to the overall transcriptomic index across developmental stages or cell
types.

## Usage

``` r
plot_contribution(phyex_set, type = c("stacked", "line"))
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet)

- type:

  One of "stacked" or "line". If stacked, will show the lines stacked in
  a cumulative area plot, otherwise they will not be stacked.

## Value

A ggplot2 object showing phylostratum contributions as a stacked area
plot

## Details

This function visualizes how different phylostrata contribute to the
overall transcriptomic index pattern across developmental stages (bulk
data) or cell types (single-cell data). Each area represents the
contribution of a specific phylostratum, with older strata typically
shown in darker colors and younger strata in lighter colors.

The plot uses the PS_colours function to create a consistent color
scheme that matches other myTAI visualizations.

## Examples

``` r
# Create contribution plot for bulk data
contrib_plot <- plot_contribution(example_phyex_set)

# Create contribution plot for single-cell data
contrib_plot_sc <- plot_contribution(example_phyex_set_sc)
```
