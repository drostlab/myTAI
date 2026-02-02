# Adaptive TXI calculation for single cell expression

Automatically selects the best TXI implementation based on dataset size
and available computational resources. Uses R implementation for smaller
datasets and C++ implementation with optimal parallelization for larger
datasets.

## Usage

``` r
.TXI_sc_adaptive(expression, strata_values, force_method = NULL, ncores = NULL)
```

## Arguments

- expression:

  Matrix of expression values, dgCMatrix

- strata_values:

  Numeric vector of phylostratum values

- force_method:

  Character, force specific method: "r", "cpp_simple", or "cpp_batched"

- ncores:

  Integer, number of cores to use (default: parallel::detectCores())

## Value

Vector of TXI values

## Details

Based on performance benchmarking: - R implementation is fastest for \<
50,000 cells - C++ batched implementation becomes advantageous for \>=
100,000 cells - Optimal core count scales with dataset size
