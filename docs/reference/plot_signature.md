# Plot Transcriptomic Signature

Create a plot of the transcriptomic index signature across developmental
stages or cell types, with options for showing individual samples/cells
and statistical testing.

## Usage

``` r
plot_signature(
  phyex_set,
  show_reps = TRUE,
  show_p_val = TRUE,
  conservation_test = stat_flatline_test,
  colour = NULL,
  ...
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet)

- show_reps:

  Logical, whether to show individual replicates

- show_p_val:

  Logical, whether to show conservation test p-value

- conservation_test:

  Function, conservation test to use for p-value calculation

- colour:

  Character, custom color for the plot elements

- ...:

  Additional arguments passed to specific methods

## Value

A ggplot2 object showing the transcriptomic signature

## Details

This function creates visualizations appropriate for the data type:

\*\*Bulk data (BulkPhyloExpressionSet):\*\* - Line plots showing TXI
trends across developmental stages - Optional individual biological
replicates as jittered points - Optional conservation test p-values

\*\*Single-cell data (ScPhyloExpressionSet):\*\* - Sina plots showing
TXI distributions across cell types or other identities - Mean TXI
values overlaid as line - Optional individual cells using geom_sina for
better visualization - Flexible identity selection from metadata via
additional parameters: - \`primary_identity\`: Character, name of
metadata column for x-axis (default: current selected identities) -
\`secondary_identity\`: Character, name of metadata column for
coloring/faceting - \`facet_by_secondary\`: Logical, whether to facet by
secondary identity (default: FALSE uses colouring)

## Examples

``` r
# Basic signature plot for bulk data
p <- plot_signature(example_phyex_set)
#> 
Computing: [========================================] 100% (done)                         
#> Ran Flat Line Test
#> Significance status of signature: not significant (= no evolutionary signature in the transcriptome).
```
