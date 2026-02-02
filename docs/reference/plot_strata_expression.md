# Plot Expression Levels by Phylostratum

Create a boxplot showing the distribution of expression levels for each
phylostratum.

## Usage

``` r
plot_strata_expression(phyex_set, aggregate_FUN = mean)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet)

- aggregate_FUN:

  Function to aggregate expression across identities (default: mean)

## Value

A ggplot2 object showing expression distributions by phylostratum

## Details

This function creates a boxplot visualization showing how expression
levels vary across different phylostrata. Each point represents a gene,
and the boxes show the distribution of expression levels within each
phylostratum.

For bulk data, expression is aggregated across developmental stages. For
single-cell data, expression is aggregated across cell types.

## Examples

``` r
# Plot expression by strata using mean aggregation for bulk data
p1 <- plot_strata_expression(example_phyex_set, aggregate_FUN = mean)

# Plot using median aggregation for single-cell data
p2 <- plot_strata_expression(example_phyex_set_sc, aggregate_FUN = median)
```
