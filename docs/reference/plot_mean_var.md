# Plot Mean-Variance Relationship

Create a scatter plot showing the relationship between mean expression
and variance for genes, colored by phylostratum, with optional
highlighting and labeling of specific genes.

## Usage

``` r
plot_mean_var(
  phyex_set,
  highlight_genes = NULL,
  colour_by = c("none", "strata")
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet) containing gene expression data.

- highlight_genes:

  Optional character vector of gene IDs to highlight and label on the
  plot.

- colour_by:

  Character string specifying coloring scheme: "none" (default),
  "strata" colors by phylostratum

## Value

A ggplot2 object showing the mean-variance relationship.

## Details

This function plots the mean expression versus variance for each gene,
with points colored by phylostratum. Optionally, specific genes can be
highlighted and labeled. This visualization helps identify expression
patterns and heteroscedasticity in the data.

The function uses collapsed expression data (averaged across replicates
for bulk data, or averaged across cells per cell type for single-cell
data).

## Examples

``` r
# Create mean-variance plot for bulk data
mv_plot <- plot_mean_var(example_phyex_set)

# Highlight and label specific genes in single-cell data
mv_plot_sc <- plot_mean_var(example_phyex_set_sc, 
  highlight_genes = example_phyex_set_sc@gene_ids[1:3])

# Color by phylostratum
mv_plot_colored <- plot_mean_var(example_phyex_set, colour_by = "strata")
```
