# Diagnose Test Robustness

Evaluate the robustness of conservation tests across different sample
sizes for null distribution generation.

## Usage

``` r
diagnose_test_robustness(
  test,
  phyex_set,
  sample_sizes = c(500, 1000, 5000, 10000),
  plot_result = TRUE,
  num_reps = 5,
  ...
)
```

## Arguments

- test:

  Function representing the conservation test to evaluate

- phyex_set:

  A PhyloExpressionSet object

- sample_sizes:

  Numeric vector of sample sizes to test (default: c(500, 1000, 5000,
  10000))

- plot_result:

  Logical indicating whether to plot results (default: TRUE)

- num_reps:

  Number of replicates for each sample size (default: 5)

- ...:

  Additional arguments passed to the test function

## Value

A data frame with test results across different sample sizes

## Details

This function assesses how consistent test results are across different
sample sizes for null distribution generation, helping to determine
appropriate sample sizes for reliable testing.

## Examples

``` r
# Diagnose flatline test robustness
p <- example_phyex_set
robustness <- diagnose_test_robustness(stat_flatline_test, 
                                       p |> select_genes(p@gene_ids[1:30]),
                                       sample_sizes=c(10,20),
                                       plot_result=FALSE,
                                       num_reps=2)
#> 
Computing: [========================================] 100% (done)                         
#> 
Computing: [========================================] 100% (done)                         
#> 
Computing: [========================================] 100% (done)                         
#> 
Computing: [========================================] 100% (done)                         
#> 
Computing: [========================================] 100% (done)                         
```
