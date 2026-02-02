# Calculate Consensus Gene Set

Calculate consensus genes from multiple GATAI runs based on frequency
threshold.

## Usage

``` r
consensus(x, p = 0.5)
```

## Arguments

- x:

  List of gene sets from different GATAI runs

- p:

  Frequency threshold (default: 0.5)

## Value

Character vector of consensus genes

## Details

This function identifies genes that appear in at least p proportion of
the input gene sets, providing a consensus set of genes across multiple
GATAI runs. \# Calculate consensus from multiple runs \# consensus_genes
\<- consensus(gatai_runs, p = 0.5)
