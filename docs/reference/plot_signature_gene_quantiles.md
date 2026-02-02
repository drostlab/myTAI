# Plot Signature Across Gene Expression Quantiles

Create a plot showing how the transcriptomic signature changes when
genes are progressively removed based on expression quantiles.

## Usage

``` r
plot_signature_gene_quantiles(
  phyex_set,
  quantiles = c(1, 0.99, 0.95, 0.9, 0.8),
  selection_FUN = genes_top_mean,
  ...
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet)

- quantiles:

  Numeric vector of quantiles to test (default: c(1.0, 0.99, 0.95, 0.90,
  0.80))

- selection_FUN:

  Function to select genes for removal (default: genes_top_mean)

- ...:

  Additional arguments passed to plot_signature_multiple

## Value

A ggplot2 object showing signatures across different quantiles

## Details

This function systematically removes genes based on expression quantiles
and shows how the transcriptomic signature changes. This is useful for
understanding the contribution of highly expressed genes to the overall
pattern and for assessing the robustness of phylotranscriptomic
patterns.

The analysis works with both bulk and single-cell data, helping to
determine whether phylotranscriptomic patterns are driven by a few
highly expressed genes or represent broad transcriptomic trends.

## Examples

``` r
# Plot signature across expression quantiles for bulk data
phyex_set <- example_phyex_set |>
    select_genes(example_phyex_set@gene_ids[1:100])
phyex_set@null_conservation_sample_size <- 500
p <- plot_signature_gene_quantiles(phyex_set, quantiles = c(0.95, 0.90))
#> 
Computing: [========================================] 100% (done)                         
```
