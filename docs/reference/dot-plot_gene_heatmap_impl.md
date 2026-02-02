# Shared Gene Heatmap Implementation

Internal helper function that contains the shared logic for creating
gene heatmaps.

## Usage

``` r
.plot_gene_heatmap_impl(
  expression_matrix,
  strata,
  gene_ids,
  num_strata,
  genes = NULL,
  top_p = NULL,
  top_k = 30,
  std = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  order_by_stage = TRUE,
  show_gene_age = TRUE,
  show_gene_ids = FALSE,
  gene_annotation = NULL,
  gene_annotation_colors = NULL,
  annotation_col = NULL,
  annotation_col_colors = NULL,
  ...
)
```

## Arguments

- expression_matrix:

  Matrix of expression values (genes x samples)

- strata:

  Factor vector of gene phylostrata

- gene_ids:

  Character vector of all gene IDs in the dataset (used for phylostratum
  mapping)

- num_strata:

  Integer number of phylostrata

- genes:

  Character vector of specific genes to plot. If NULL, uses top dynamic
  genes

- top_p:

  Proportion of most dynamic genes to include (default: NULL). Ignored
  if top_k is specified.

- top_k:

  Absolute number of top genes to select (default: 30). Takes precedence
  over top_p.

- std:

  Logical indicating whether to use standardized expression values
  (default: TRUE)

- cluster_rows:

  Logical indicating whether to cluster genes/rows (default: FALSE)

- cluster_cols:

  Logical indicating whether to cluster identities/columns (default:
  FALSE)

- order_by_stage:

  Logical indicating whether to order genes by expression
  angle/developmental stage (default: TRUE). Ignored if cluster_rows is
  TRUE

- show_gene_age:

  Logical indicating whether to show gene age as row annotation
  (default: TRUE)

- show_gene_ids:

  Logical indicating whether to show gene names (default: FALSE)

- gene_annotation:

  Data frame with custom gene annotations, rownames should match gene
  IDs (default: NULL)

- gene_annotation_colors:

  Named list of color vectors for custom gene annotations (default:
  NULL)

- annotation_col:

  Data frame with column annotations (default: NULL)

- annotation_col_colors:

  List of colors for column annotations (default: NULL)

- ...:

  Additional arguments passed to pheatmap::pheatmap

## Value

A ggplot object showing the gene expression heatmap
