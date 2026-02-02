# Match Expression Matrix with Phylostratum Map

Join single-cell gene expression matrix with a phylostratum mapping to
create a ScPhyloExpressionSet object.

## Usage

``` r
match_map_sc_matrix(
  expression_matrix,
  metadata,
  phylomap,
  strata_legend = NULL,
  ...
)
```

## Arguments

- expression_matrix:

  Expression matrix with genes as rows and cells as columns

- metadata:

  Data frame with cell metadata, rownames should match colnames of
  expression_matrix

- phylomap:

  A data frame with two columns: phylostratum assignments and gene IDs

- strata_legend:

  A data frame with two columns: phylostratum assignments and name of
  each stratum. If NULL, numeric labels will be used (default: NULL)

- ...:

  Additional arguments passed to ScPhyloExpressionSet_from_matrix

## Value

A ScPhyloExpressionSet object

## Details

This function combines phylostratum mapping with expression matrix and
metadata to create a ScPhyloExpressionSet. Only genes present in both
the expression matrix and phylomap will be retained. All discrete
metadata columns are converted to factors automatically.
