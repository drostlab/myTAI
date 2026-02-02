# Fit Gamma Distribution Parameters

Fit gamma distribution parameters using a robust method that filters
outliers iteratively to find the best fit.

## Usage

``` r
.fit_gamma(x)
```

## Arguments

- x:

  Numeric vector of data values

## Value

List with shape and rate parameters

## Details

This function uses an iterative approach to filter outliers and find the
gamma distribution parameters that best fit the data, improving
robustness compared to standard fitting methods. \# Fit gamma
distribution \# params \<- .fit_gamma(data_vector)
