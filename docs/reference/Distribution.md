# Distribution S7 Class

S7 class for representing probability distributions used in statistical
testing, including PDF, CDF, quantile functions, and fitting procedures.

## Usage

``` r
Distribution(
  name = stop("@name is required"),
  pdf = stop("@pdf is required"),
  cdf = stop("@cdf is required"),
  quantile_function = stop("@quantile_function is required"),
  fitting_function = stop("@fitting_function is required"),
  param_names = stop("@param_names is required")
)
```

## Arguments

- name:

  Character string identifying the distribution

- pdf:

  Function for probability density function

- cdf:

  Function for cumulative distribution function

- quantile_function:

  Function for quantile calculations

- fitting_function:

  Function to fit distribution parameters from data

- param_names:

  Character vector of parameter names

## Details

The Distribution class provides a unified interface for different
probability distributions used in phylotranscriptomic testing. Each
distribution includes the necessary functions for statistical inference.

## Examples

``` r
# Access predefined distributions
normal_dist <- distributions$normal
gamma_dist <- distributions$gamma
```
