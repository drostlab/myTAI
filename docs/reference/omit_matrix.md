# Compute TXI Profiles Omitting Each Gene

For each gene i, exclude the corresponding gene i from the
PhyloExpressionSet and compute the TXI profile for the dataset with gene
i excluded.

This procedure results in a TXI profile matrix storing the TXI profile
for each omitted gene i.

## Usage

``` r
omit_matrix(phyex_set)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

## Value

A numeric matrix storing TXI profiles for each omitted gene i

## Details

This function systematically removes each gene and recalculates the
transcriptomic index profile to assess the contribution of individual
genes to the overall pattern. This is useful for identifying genes that
have a large influence on the phylotranscriptomic signature.

## Author

Hajk-Georg Drost

## Examples

``` r
# Compute omit matrix for a PhyloExpressionSet
omit_mat <- omit_matrix(example_phyex_set)
```
