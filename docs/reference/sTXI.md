# Calculate Stratum-Specific Transcriptomic Index

Calculate the stratum-specific transcriptomic index (sTXI) by summing
pTXI values within each phylostratum.

## Usage

``` r
sTXI(phyex_set, option = "identity")
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- option:

  Character string specifying calculation method: - "identity": Sum pTXI
  values within each stratum - "add": Cumulative sum across strata

## Value

Matrix of sTXI values with strata as rows and identities as columns

## Examples

``` r
# Calculate sTXI values
stxi_values <- sTXI(example_phyex_set, option = "identity")
stxi_cumsum <- sTXI(example_phyex_set, option = "add")
```
