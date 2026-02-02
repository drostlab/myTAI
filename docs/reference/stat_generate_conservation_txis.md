# Generate Null Conservation TXI Distribution

Generate a null distribution of transcriptomic index values for
conservation testing by permuting phylostratum assignments.

## Usage

``` r
stat_generate_conservation_txis(strata_vector, count_matrix, sample_size)
```

## Arguments

- strata_vector:

  Numeric vector of phylostratum assignments

- count_matrix:

  Matrix of expression counts

- sample_size:

  Number of permutations to generate

## Value

Matrix of permuted TXI values for null hypothesis testing

## Details

This function creates a null distribution by randomly permuting
phylostratum assignments while keeping the expression matrix fixed. This
preserves the expression structure while breaking the
phylostratum-expression relationship. \# Generate null TXI distribution
\# null_txis \<- stat_generate_conservation_txis(strata_vec,
expr_matrix, 1000)
