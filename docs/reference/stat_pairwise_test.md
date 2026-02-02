# Pairwise Conservation Test

Test for significant differences in transcriptomic index values between
two groups of developmental stages.

## Usage

``` r
stat_pairwise_test(phyex_set, modules, alternative = c("greater", "less"), ...)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- modules:

  A named list with elements 'contrast1' and 'contrast2' containing
  stage indices for each contrast group

- alternative:

  Character string specifying the alternative hypothesis: "greater"
  (contrast1 \> contrast2) or "less" (contrast1 \< contrast2)

- ...:

  Additional arguments passed to stat_generic_conservation_test

## Value

A ConservationTestResult object with pairwise test results

## Details

The pairwise test compares the mean transcriptomic index values between
two groups of developmental stages. This is useful for testing specific
hypotheses about differences in gene age composition between
developmental periods.

## See also

[`stat_generic_conservation_test`](https://drostlab.github.io/myTAI/reference/stat_generic_conservation_test.md)

## Author

Jaruwatana Sodai Lotharukpong

## Examples

``` r
# Define contrast groups
modules <- list(contrast1 = 1:3, contrast2 = 4:7)
result <- stat_pairwise_test(example_phyex_set, modules, alternative = "greater")

```
