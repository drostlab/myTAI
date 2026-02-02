# Memoized Null Conservation TXI Generation

Memoized version of generate_conservation_txis for improved performance
with repeated calls using the same parameters.

## Usage

``` r
.memo_generate_conservation_txis(strata_vector, count_matrix, sample_size)
```

## Details

This function caches results to avoid recomputing expensive permutations
when the same parameters are used multiple times.
