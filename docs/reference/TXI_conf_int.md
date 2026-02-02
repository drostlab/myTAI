# Confidence Intervals for Transcriptomic Index (TXI)

Compute confidence intervals for the TXI using bootstrapped TXI values.

## Usage

``` r
TXI_conf_int(phyex_set, probs = c(0.025, 0.975))
```

## Arguments

- phyex_set:

  A BulkPhyloExpressionSet object

- probs:

  Numeric vector of probabilities for the confidence interval (default:
  c(0.025, 0.975))

## Value

A tibble with first column Identity names, second column lower bound,
third column upper bound

## Details

This function returns confidence intervals for the TXI for each identity
(sample or group), based on the bootstrapped TXI values stored in the
PhyloExpressionSet object.
