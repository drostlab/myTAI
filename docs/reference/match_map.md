# Match Gene Expression Data with Phylostratum Map

Join gene expression data with a phylostratum mapping to create a
BulkPhyloExpressionSet object.

## Usage

``` r
match_map(
  data,
  phylomap,
  groups = colnames(data[, 2:ncol(data)]),
  name = NULL,
  ...
)
```

## Arguments

- data:

  A data frame where column 1 contains gene IDs and columns 2+ contain
  expression data

- phylomap:

  A data frame with two columns: phylostratum assignments and gene IDs

- groups:

  A factor or character vector indicating which group each sample
  belongs to. Default uses column names from expression data

- name:

  A character string naming the dataset. Default uses the variable name

- ...:

  Additional arguments passed to as_BulkPhyloExpressionSet

## Value

A BulkPhyloExpressionSet object

## Examples

``` r
# Match expression data with phylostratum map
# bulk_set <- match_map(expression_data, phylo_map, 
#                       groups = c("stage1", "stage2", "stage3"),
#                       name = "Matched Dataset")
```
