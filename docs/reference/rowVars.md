# Row-wise Variance Calculation

Calculate the variance for each row of a matrix or data frame.

## Usage

``` r
rowVars(x, ...)
```

## Arguments

- x:

  A numeric matrix or data frame

- ...:

  Additional arguments passed to rowSums and rowMeans

## Value

A numeric vector containing the variance for each row

## Details

This function computes the sample variance for each row using the
formula: var = sum((x - mean(x))^2) / (n - 1) \# Calculate row variances
for a matrix \# mat \<- matrix(1:12, nrow = 3) \# row_vars \<-
rowVars(mat)
