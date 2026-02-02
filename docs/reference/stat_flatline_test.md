# Flat Line Test for Conservation Pattern

Perform a flat line test to assess whether the transcriptomic index
profile shows a flat (non-varying) pattern across developmental stages.

## Usage

``` r
stat_flatline_test(phyex_set, ...)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- ...:

  Additional arguments passed to stat_generic_conservation_test

## Value

A test result object containing p-value and test statistics

## Details

The flat line test evaluates whether the TXI profile remains constant
across developmental stages by testing the variance of the profile
against a null distribution. A significant result indicates rejection of
the flat line pattern.

## See also

[`stat_generic_conservation_test`](https://drostlab.github.io/myTAI/reference/stat_generic_conservation_test.md)

## Examples

``` r
# Perform flat line test
result <- stat_flatline_test(example_phyex_set)

```
