# Match Single-Cell Expression Data with Phylostratum Map (Seurat)

Join single-cell gene expression data (from a Seurat object) with a
phylostratum mapping to create a ScPhyloExpressionSet object.
Automatically extracts dimensional reductions and metadata.

## Usage

``` r
match_map_sc_seurat(
  seurat,
  phylomap,
  layer = "counts",
  strata_legend = NULL,
  ...
)
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
  each stratum. If NULL, numeric labels will be used (default: NULL)

- ...:

  Additional arguments passed to ScPhyloExpressionSet_from_seurat

## Value

A ScPhyloExpressionSet object

## Details

This is a convenience function that combines phylostratum mapping with
Seurat object conversion. Only genes present in both the expression data
and phylomap will be retained. The function extracts all metadata and
dimensional reductions from the Seurat object.
