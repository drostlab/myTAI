# Test Result S7 Class

S7 class for storing and manipulating statistical test results from
phylotranscriptomic conservation tests.

## Usage

``` r
TestResult(
  method_name = stop("@method_name is required"),
  test_stat = stop("@test_stat is required"),
  fitting_dist = stop("@fitting_dist is required"),
  params = stop("@params is required"),
  alternative = "two-sided",
  null_sample = stop("@null_sample is required"),
  data_name = character(0),
  p_label = "p_val"
)
```

## Arguments

- method_name:

  Character string identifying the test method

- test_stat:

  Numeric test statistic value

- fitting_dist:

  Distribution object used for null hypothesis testing

- params:

  List of fitted distribution parameters

- alternative:

  Character string specifying alternative hypothesis ("two-sided",
  "less", "greater")

- null_sample:

  Numeric vector of null distribution samples

- data_name:

  Character string naming the dataset (optional)

- p_label:

  Character string for p-value label (default: "p_val")

## Value

A TestResult object

## Details

The TestResult class provides computed properties including: -
\`p_value\`: Computed p-value based on test statistic and fitted
distribution
