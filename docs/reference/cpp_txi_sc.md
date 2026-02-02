# Calculate TXI for Single-Cell Expression Data (C++ Implementation)

Efficiently calculate TXI values for sparse single-cell expression
matrices using batch processing and parallel computation.

## Usage

``` r
cpp_txi_sc(expression, strata_values, batch_size = 2000L, ncores = 10L)
```

## Arguments

- expression:

  Sparse expression matrix (genes x cells) - dgCMatrix format

- strata_values:

  Numeric vector of phylostratum values for each gene

- batch_size:

  Integer, number of cells to process per batch (default: 2000)

- ncores:

  Integer, number of cores to use for parallel processing (default: 10,
  automatically capped at available cores)

## Value

Numeric vector of TXI values for each cell

## Details

This function processes large sparse single-cell expression matrices
efficiently by: - Splitting cells into batches to manage memory usage -
Using parallel processing across batches - Leveraging sparse matrix
operations to skip zero entries - Handling cells with zero expression by
returning NA

## Author

Kristian K Ullrich
