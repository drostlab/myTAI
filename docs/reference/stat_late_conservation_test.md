# Late Conservation Test

Test for late conservation patterns in transcriptomic data by comparing
late developmental stages to early and mid stages.

## Usage

``` r
stat_late_conservation_test(phyex_set, modules, ...)
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

A ConservationTestResult object with late conservation test results

## Details

The late conservation test evaluates whether later developmental stages
show lower transcriptomic index values (indicating older genes) compared
to earlier stages. The test computes a score based on the minimum
difference between early vs. late and mid vs. late TXI values.

## See also

[`stat_generic_conservation_test`](https://drostlab.github.io/myTAI/reference/stat_generic_conservation_test.md),
[`stat_early_conservation_test`](https://drostlab.github.io/myTAI/reference/stat_early_conservation_test.md)

## Examples

``` r
# Define developmental modules
modules <- list(early = 1:2, mid = 3:5, late = 6:7)
result <- stat_late_conservation_test(example_phyex_set_old, modules)
#> 
Computing: [========================================] 100% (done)                         

```
