# Plot Distribution of Genes Across Phylostrata

Create a bar plot showing the distribution of genes across phylostrata,
with options for showing observed vs. expected ratios.

## Usage

``` r
plot_distribution_strata(
  strata,
  selected_gene_ids = names(strata),
  as_log_obs_exp = FALSE
)
```

## Arguments

- strata:

  Named factor vector of phylostratum assignments (names are gene IDs)

- selected_gene_ids:

  Character vector of gene IDs to include in the plot (default: all
  genes in strata)

- as_log_obs_exp:

  Logical indicating whether to show log2(observed/expected) ratios
  instead of raw counts (default: FALSE)

## Value

A ggplot2 object showing the phylostratum distribution

## Details

This function visualizes how genes are distributed across different
phylostrata. When as_log_obs_exp=FALSE, it shows raw gene counts per
stratum. When TRUE, it shows log2 ratios of observed vs. expected gene
counts, useful for identifying enrichment or depletion of specific
strata in gene sets.

## Examples

``` r
# Plot raw gene counts by strata
p1 <- plot_distribution_strata(example_phyex_set@strata)

# Plot observed vs expected ratios for selected genes
p2 <- plot_distribution_strata(example_phyex_set@strata, 
                               selected_gene_ids = example_phyex_set@gene_ids[5:20],
                               as_log_obs_exp = TRUE)
#> Warning: Ignoring unknown parameters: `size`
```
