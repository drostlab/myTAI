# Plot Individual Gene Expression Profiles

Create a plot showing expression profiles for individual genes across
developmental stages or cell types, with various visualization options.

## Usage

``` r
plot_gene_profiles(
  phyex_set,
  genes = NULL,
  show_set_mean = FALSE,
  show_reps = FALSE,
  transformation = c("log", "std_log", "none"),
  colour_by = c("manual", "strata", "stage"),
  colours = NULL,
  max_genes = 100,
  show_labels = TRUE,
  label_size = 1.75,
  show_legend = TRUE,
  facet_by_strata = FALSE
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet)

- genes:

  Character vector of gene IDs to plot. If NULL, top expressing genes
  are selected

- show_set_mean:

  Logical indicating whether to show the mean expression across all
  genes (default: FALSE)

- show_reps:

  Logical indicating whether to show individual replicates (bulk) or
  cells (single-cell) (default: FALSE)

- transformation:

  Character string specifying expression transformation: "log" (log1p),
  "std_log" (standardized log1p), or "none" (default: "log")

- colour_by:

  Character string specifying coloring scheme: "strata" (by
  phylostratum), "stage" (by developmental stage/cell type), or "manual"
  (default: "manual")

- colours:

  Optional vector of colors for manual coloring (default: NULL)

- max_genes:

  Maximum number of genes to plot when genes=NULL (default: 100)

- show_labels:

  Logical indicating whether to show gene labels (default: TRUE)

- label_size:

  Font size of gene id labels if shown (default: 0.5).

- show_legend:

  Logical indicating whether to show legend (default: TRUE)

- facet_by_strata:

  Logical indicating whether to facet by phylostratum (default: FALSE)

## Value

A ggplot2 object showing gene expression profiles

## Details

This function creates detailed visualizations of individual gene
expression patterns across development (bulk data) or cell types
(single-cell data). Genes can be colored by phylostratum or
developmental stage, and various transformations can be applied to
highlight different aspects of the data.

## Author

Filipa Martins Costa, Stefan Manolache, Hajk-Georg Drost

## Examples

``` r
# Plot specific genes for bulk data
p1 <- plot_gene_profiles(example_phyex_set, genes = example_phyex_set@gene_ids[1:5])

# Plot for single-cell data with faceting by strata
p2 <- plot_gene_profiles(example_phyex_set_sc, facet_by_strata = TRUE)
```
