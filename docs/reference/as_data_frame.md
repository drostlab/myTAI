# Convert BulkPhyloExpressionSet to Data Frame

Convert a BulkPhyloExpressionSet object back to the original data frame
format with phylostratum, gene ID, and expression data as columns.

## Usage

``` r
as_data_frame(phyex_set, use_collapsed = FALSE)
```

## Arguments

- phyex_set:

  A BulkPhyloExpressionSet object

- use_collapsed:

  Logical indicating whether to use collapsed expression data (default:
  FALSE)

## Value

A data frame where column 1 contains phylostratum information, column 2
contains gene IDs, and columns 3+ contain expression data

## Examples

``` r
# Convert BulkPhyloExpressionSet back to data frame
df <- as_data_frame(example_phyex_set)
df_collapsed <- as_data_frame(example_phyex_set, use_collapsed = TRUE)
```
