# Downsample ScPhyloExpressionSet

Create a downsampled copy of a ScPhyloExpressionSet object with fewer
cells per identity.

## Usage

``` r
downsample(phyex_set, downsample = 10)
```

## Arguments

- phyex_set:

  A ScPhyloExpressionSet object

- downsample:

  Integer, number of cells to keep per identity (default: 10)

## Value

A new ScPhyloExpressionSet object with downsampled cells

## Details

This function creates a new ScPhyloExpressionSet with a subset of cells,
maintaining the same proportional representation across identities. The
sampling is stratified by the current `selected_idents` grouping. All
metadata and reductions are filtered to match the selected cells.

## Examples

``` r
# Downsample to 20 cells per identity
small_set <- downsample(example_phyex_set_sc, downsample = 20)
#> Downsampled from 200 to 60 cells.

# Change grouping and downsample
example_phyex_set_sc@selected_idents <- "day"
treatment_set <- downsample(example_phyex_set_sc, downsample = 15)
#> Downsampled from 200 to 60 cells.
```
