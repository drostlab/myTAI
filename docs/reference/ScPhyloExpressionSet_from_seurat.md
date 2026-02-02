# Convert Seurat Object to Single-Cell PhyloExpressionSet

Convert a Seurat object with phylostratum information into a
ScPhyloExpressionSet object for single-cell phylotranscriptomic
analysis. Automatically extracts dimensional reductions if present, or
computes basic PCA and UMAP if none are available.

## Usage

``` r
ScPhyloExpressionSet_from_seurat(
  seurat,
  strata,
  layer = "counts",
  selected_idents = NULL,
  name = "Single-Cell PhyloExpressionSet",
  seed = 42,
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

- selected_idents:

  Character string specifying which metadata column to use for grouping
  (default: NULL, uses active idents)

- name:

  A character string naming the dataset (default: "Single-Cell Phylo
  Expression Set")

- seed:

  Integer seed for reproducible UMAP computation (default: 42)

- ...:

  Additional arguments passed to ScPhyloExpressionSet constructor

## Value

A ScPhyloExpressionSet object
