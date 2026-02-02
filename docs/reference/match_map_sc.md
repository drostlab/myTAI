# Match Single-Cell Expression Data with Phylostratum Map

Join single-cell gene expression data (from a Seurat object) with a
phylostratum mapping to create a ScPhyloExpressionSet object.

## Usage

``` r
match_map_sc(seurat, phylomap, layer = "counts", strata_legend = NULL, ...)
```

## Arguments

- seurat:

  A Seurat object containing single-cell expression data

- phylomap:

  A data frame with two columns: phylostratum assignments and gene IDs

- layer:

  Character string specifying which layer to use from the Seurat object
  (default: "counts")

- strata_legend:

  A data frame with two columns: phylostratum assignments and name of
  each stratum. If NULL, no labels will be added (default: NULL) If
  NULL, uses sorted unique values from column 1

- ...:

  Additional arguments passed to as_ScPhyloExpressionSet

## Value

A ScPhyloExpressionSet object

## Examples

``` r
# Match Seurat object with phylostratum map
# sc_set <- match_map_sc(seurat_obj, phylo_map, layer = "counts", name = "SC Matched Dataset")
```
