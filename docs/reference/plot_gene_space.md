# Plot Gene Space Using PCA

Create a PCA plot showing genes in expression space with ideal
expression patterns overlaid as reference points.

## Usage

``` r
plot_gene_space(
  phyex_set,
  top_p = 0.2,
  genes = NULL,
  colour_by = c("identity", "strata")
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet)

- top_p:

  Proportion of most dynamic genes to include when genes=NULL (default:
  0.2)

- genes:

  Character vector of specific genes to plot. If NULL, uses top dynamic
  genes

- colour_by:

  Character string specifying coloring scheme: "identity" (by peak
  expression stage/identity) or "strata" (by phylostratum) (default:
  "identity")

## Value

A ggplot2 object showing the gene space PCA plot

## Details

This function creates a PCA visualization of genes in expression space,
with ideal expression patterns (early, mid, late, reverse mid) overlaid
as reference points. The analysis uses log-transformed and standardized
expression values. Genes are colored either by their phylostratum or by
their peak expression stage.

## Examples

``` r
# Plot gene space colored by identity
p1 <- plot_gene_space(example_phyex_set, colour_by = "identity")

# Plot specific genes colored by strata
p2 <- plot_gene_space(example_phyex_set, 
                      genes = example_phyex_set@gene_ids[1:5], 
                      colour_by = "strata")
```
