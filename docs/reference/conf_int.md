# Calculate Confidence Intervals for Test Result

Calculate confidence intervals for the null distribution of a test
result.

## Usage

``` r
conf_int(test_result, probs = c(0.025, 0.975))
```

## Arguments

- test_result:

  A TestResult object

- probs:

  Numeric vector of probabilities for quantiles (default: c(0.025,
  0.975))

## Value

Numeric vector of quantiles from the null distribution

## Details

\# Calculate 95 \# ci \<- conf_int(test_result, probs = c(0.025, 0.975))
