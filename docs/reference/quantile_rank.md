# Calculate Quantile Ranks

Calculate quantile ranks for a numeric vector, handling ties using
average method.

## Usage

``` r
quantile_rank(x)
```

## Arguments

- x:

  numeric vector for which to calculate quantile ranks

## Value

A numeric vector of quantile ranks between 0 and 1

## Examples

``` r
# Calculate quantile ranks for a vector
ranks <- quantile_rank(c(1, 2, 3, 4, 5))
```
