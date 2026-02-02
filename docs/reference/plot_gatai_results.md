# Plot Comprehensive GATAI Results

Create a suite of plots summarizing the effects of GATAI gene removal on
phylotranscriptomic patterns.

## Usage

``` r
plot_gatai_results(
  phyex_set,
  gatai_result,
  conservation_test = stat_flatline_test,
  runs_threshold = 0.5,
  signature_plot_type = c("separate", "combined")
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object containing the original gene expression
  data.

- gatai_result:

  Result list from
  [`destroy_pattern()`](https://drostlab.github.io/myTAI/reference/destroy_pattern.md),
  containing GATAI analysis output.

- conservation_test:

  Function for conservation test (default: `stat_flatline_test`).

- runs_threshold:

  Threshold for gene removal consistency across runs (default: 0.5).

- signature_plot_type:

  Type of signature plot: "separate" for individual plots, "combined"
  for overlay (default: both options).

## Value

A named list of ggplot/patchwork objects and results:

- signature_plots:

  Signature plots before/after GATAI and top variance removal

- heatmap_plot:

  Heatmap of GATAI-removed genes

- profiles_plot:

  Gene expression profiles of GATAI-removed genes

- profiles_plot_facet:

  Faceted gene profiles by strata

- gene_space_plot:

  Gene space plot of GATAI-removed genes

- mean_var_plot:

  Mean-variance plot highlighting GATAI-removed genes

- strata_plot:

  Phylostrata distribution plot (log obs/exp) for GATAI-removed genes

- null_dist_plot:

  Null distribution plot with test statistics and p-values

- convergence_plots:

  GATAI convergence plots (if available)

## Details

This function provides a comprehensive visualization of the impact of
GATAI gene removal, including transcriptomic signature plots, gene
expression profiles, heatmaps, mean-variance relationships, phylostrata
distributions, conservation test comparisons, and convergence
diagnostics.

## Author

Filipa Martins Costa, Stefan Manolache
