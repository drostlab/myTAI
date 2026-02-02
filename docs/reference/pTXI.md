# Calculate Phylostratum-Specific Transcriptomic Index

Calculate pTXI values for expression data. This is a generic function
that dispatches based on the input type.

## Usage

``` r
pTXI(phyex_set, reps = FALSE)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- reps:

  Whether to return pTXI for each sample instead for each group

## Value

Matrix of pTXI values
