# Get Strata Legend from PhyloExpressionSet

Extract the strata legend (phylostratum ranks and names) from a
PhyloExpressionSet object.

## Usage

``` r
get_strata_legend(phyex_set)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

## Value

A tibble with two columns: Rank (phylostratum values) and Name
(phylostratum labels), sorted by Rank

## Examples

``` r
# Get strata legend
legend <- get_strata_legend(example_phyex_set)
```
