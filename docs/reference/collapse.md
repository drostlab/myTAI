# Collapse PhyloExpressionSet Replicates

Convert a PhyloExpressionSet with replicates to one with collapsed
expression data.

## Usage

``` r
collapse(phyex_set, ...)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- ...:

  Additional arguments passed to methods

## Value

A new PhyloExpressionSet object with collapsed expression data

## Examples

``` r
# Collapse replicates in a PhyloExpressionSet
collapsed_set <- collapse(example_phyex_set)
```
