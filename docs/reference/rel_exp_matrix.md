# Compute Relative Expression Matrix for PhyloExpressionSet

Computes relative expression profiles for all age categories in a
PhyloExpressionSet.

## Usage

``` r
rel_exp_matrix(phyex_set)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object (BulkPhyloExpressionSet or
  ScPhyloExpressionSet).

## Value

A matrix with age categories as rows and identities as columns,
containing relative expression values.
