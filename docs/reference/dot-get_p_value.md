# Calculate P-Value from Distribution

Internal function to calculate p-values from cumulative distribution
functions based on test statistics and alternative hypothesis
specifications.

## Usage

``` r
.get_p_value(
  cdf,
  test_stat,
  params,
  alternative = c("two-sided", "less", "greater")
)
```

## Arguments

- cdf:

  Cumulative distribution function

- test_stat:

  Numeric test statistic value

- params:

  List of distribution parameters

- alternative:

  Character string specifying alternative hypothesis ("two-sided",
  "less", "greater")

## Value

Numeric p-value

## Details

This function calculates p-values using the appropriate tail(s) of the
distribution: - "greater": Uses upper tail (1 - CDF) - "less": Uses
lower tail (CDF) - "two-sided": Uses 2 \* minimum of both tails \#
Calculate p-value (internal use) \# pval \<- .get_p_value(pnorm, 1.96,
list(mean=0, sd=1), "two-sided")
