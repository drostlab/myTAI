# Convert Seurat Object to Single-Cell PhyloExpressionSet

Convert a Seurat object with phylostratum information into a
ScPhyloExpressionSet object for single-cell phylotranscriptomic
analysis.

## Usage

``` r
as_ScPhyloExpressionSet(
  seurat,
  strata,
  layer = "counts",
  name = "Single-Cell Phylo Expression Set",
  ...
)
```

## Arguments

- seurat:

  A Seurat object containing single-cell expression data

- strata:

  Factor vector of phylostratum assignments for each gene

- layer:

  Character string specifying which layer to use from the Seurat object
  (default: "counts")

- name:

  A character string naming the dataset (default: "Single-Cell Phylo
  Expression Set")

- ...:

  Additional arguments passed to ScPhyloExpressionSet constructor

## Value

A ScPhyloExpressionSet object

## Examples

``` r
# Convert Seurat object to ScPhyloExpressionSet
# sc_set <- as_ScPhyloExpressionSet(seurat_obj, phylo_map, "celltype",
#                                  name = "Brain Development SC Dataset")
```
