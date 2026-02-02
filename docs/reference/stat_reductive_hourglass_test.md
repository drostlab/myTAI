# Reductive Hourglass Test

Test for reductive hourglass patterns in transcriptomic data by
comparing early and late developmental stages to mid developmental
stages.

## Usage

``` r
stat_reductive_hourglass_test(phyex_set, modules, ...)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- modules:

  A named list with elements 'early', 'mid', and 'late' containing stage
  indices for each developmental module

- ...:

  Additional arguments passed to stat_generic_conservation_test

## Value

A ConservationTestResult object with reductive hourglass test results

## Details

The reductive hourglass test evaluates whether mid developmental stages
show lower transcriptomic index values (indicating older genes) compared
to both early and late stages. This creates an hourglass-shaped pattern
where ancient genes dominate during mid-development. The test computes a
score based on the minimum difference between early vs. mid and late vs.
mid TXI values.

## See also

[`stat_generic_conservation_test`](https://drostlab.github.io/myTAI/reference/stat_generic_conservation_test.md),
[`stat_reverse_hourglass_test`](https://drostlab.github.io/myTAI/reference/stat_reverse_hourglass_test.md)

## Examples

``` r
# Define developmental modules
modules <- list(early = 1:2, mid = 3:5, late = 6:7)
result <- stat_reductive_hourglass_test(example_phyex_set_old, modules=modules)

```
