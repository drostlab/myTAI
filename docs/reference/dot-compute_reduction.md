# Compute Dimensional Reduction

Compute PCA or UMAP on expression data when not available in stored
reductions.

## Usage

``` r
.compute_reduction(expression, method = c("PCA", "UMAP"), seed = 42)
```

## Arguments

- expression:

  Expression matrix with genes as rows and cells as columns

- method:

  Character string: "PCA" or "UMAP"

- seed:

  Integer seed for reproducible results (default: 42)

## Value

Matrix with cells as rows and dimensions as columns
