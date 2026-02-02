# Goodness of Fit Test

Perform a Kolmogorov-Smirnov test to assess goodness of fit between the
null sample and fitted distribution.

## Usage

``` r
goodness_of_fit(test_result)
```

## Arguments

- test_result:

  A TestResult object

## Value

A ks.test result object

## Details

This function tests whether the null sample follows the fitted
distribution using the Kolmogorov-Smirnov test. A significant result
indicates poor fit. \# Test goodness of fit \# gof_result \<-
goodness_of_fit(test_result)
