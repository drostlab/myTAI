# Normalise Stage Expression Data

Normalise expression data to a specified total expression level per
sample.

## Usage

``` r
normalise_stage_expression(phyex_set, total = 1e+06)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- total:

  Numeric value to normalise each sample to (default: 1e6)

## Value

A PhyloExpressionSet object with normalised expression data

## Examples

``` r
# Normalise to 1 million total expression per sample
normalised_set <- normalise_stage_expression(example_phyex_set, total = 1e6)
```
