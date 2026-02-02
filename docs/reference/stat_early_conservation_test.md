# Early Conservation Test

Test for early conservation patterns in transcriptomic data by comparing
early developmental stages to mid and late stages.

## Usage

``` r
stat_early_conservation_test(phyex_set, modules, ...)
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

A ConservationTestResult object with early conservation test results

## Details

The early conservation test evaluates whether early developmental stages
show lower transcriptomic index values (indicating older genes) compared
to later stages. The test computes a score based on the minimum
difference between mid vs. early and late vs. early TXI values.

## See also

[`stat_generic_conservation_test`](https://drostlab.github.io/myTAI/reference/stat_generic_conservation_test.md),
[`stat_late_conservation_test`](https://drostlab.github.io/myTAI/reference/stat_late_conservation_test.md)

## Examples

``` r
# Define developmental modules
p <- example_phyex_set_old |> 
     select_genes(example_phyex_set_old@gene_ids[1:100])
modules <- list(early = 1:2, mid = 3:5, late = 6:7)
result <- stat_early_conservation_test(p, modules=modules)
#> 
Computing: [========================================] 100% (done)                         

```
