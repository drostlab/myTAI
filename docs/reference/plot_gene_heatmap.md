# Plot Gene Expression Heatmap

Create a heatmap showing gene expression patterns across conditions with
optional dendrograms and gene age annotation.

## Usage

``` r
plot_gene_heatmap(
  phyex_set,
  genes = NULL,
  top_p = NULL,
  top_k = 30,
  std = TRUE,
  show_reps = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  order_by_stage = TRUE,
  show_gene_age = TRUE,
  show_gene_ids = FALSE,
  gene_annotation = NULL,
  gene_annotation_colors = NULL,
  ...
)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet)

- genes:

  Character vector of specific gene IDs to include in the heatmap
  (default: NULL for auto-selection of dynamic genes)

- top_p:

  Numeric value specifying the top proportion of genes to include
  (default: NULL). Ignored if top_k is specified.

- top_k:

  Absolute number of top genes to select (default: 30). Takes precedence
  over top_p.

- std:

  Logical indicating whether to standardize expression values (default:
  TRUE)

- show_reps:

  Logical indicating whether to show replicates or collapsed data
  (default: FALSE)

- cluster_rows:

  Logical indicating whether to cluster genes/rows (default: FALSE)

- cluster_cols:

  Logical indicating whether to cluster conditions/columns (default:
  FALSE)

- order_by_stage:

  Logical indicating whether to order genes by expression
  angle/developmental stage (default: TRUE). Ignored if cluster_rows is
  TRUE

- show_gene_age:

  Logical indicating whether to show gene age annotation (default: TRUE)

- show_gene_ids:

  Logical indicating whether to show gene identifiers (default: FALSE)

- gene_annotation:

  Data frame with custom gene annotations, rownames should match gene
  IDs (default: NULL)

- gene_annotation_colors:

  Named list of color vectors for custom gene annotations (default:
  NULL)

- ...:

  Additional arguments passed to specific methods

## Value

A ggplot object (converted from pheatmap) showing the gene expression
heatmap

## Details

This function creates a comprehensive heatmap visualization of gene
expression patterns. By default, genes are ordered by their expression
angle (developmental trajectory). The function supports clustering of
both genes and identities, and can optionally display gene age
(phylostratum) as a colored annotation bar.

For bulk data, the heatmap shows expression across developmental
conditions. For single-cell data, the heatmap shows expression across
cell types.

The gene age annotation uses the PS_colours function to create a
consistent color scheme across different myTAI visualizations.

Custom gene annotations can be provided via the `gene_annotation`
parameter, which should be a data frame with gene IDs as rownames and
annotation categories as columns. Corresponding colors should be
provided via `gene_annotation_colors` as a named list where names match
the annotation column names.

## Examples

``` r
# Basic heatmap with gene age annotation
p1 <- plot_gene_heatmap(example_phyex_set, show_gene_age = TRUE)

# Single-cell heatmap with subset of cells
p2 <- plot_gene_heatmap(example_phyex_set_sc, show_reps = TRUE, max_cells_per_type = 3)
#> Showing 9 out of 200 cells. Use max_cells_per_type to control the number of cells per type.

# Custom gene annotation example
gene_ids <- example_phyex_set@gene_ids[1:3]
gene_annot <- data.frame(
  Category = c("High", "Medium", "Low"),
  row.names = gene_ids
)
colors <- list(Category = c("High" = "red", "Medium" = "yellow", "Low" = "blue"))
p3 <- plot_gene_heatmap(example_phyex_set |> select_genes(gene_ids), gene_annotation = gene_annot, 
                       gene_annotation_colors = colors, show_gene_age = FALSE)
```
