# Create Full GATAI Convergence Plot

Create a comprehensive plot showing GATAI convergence across multiple
runs and thresholds.

## Usage

``` r
full_gatai_convergence_plot(phyex_set, runs, p = 0.5, ps = c(0.25, 0.5, 0.75))
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- runs:

  List of GATAI run results

- p:

  Consensus threshold for petal plot (default: 0.5)

- ps:

  Vector of consensus thresholds for convergence plots (default: c(0.25,
  0.5, 0.75))

## Value

A patchwork composition showing convergence analysis

## Details

This function creates a comprehensive visualization of GATAI convergence
including consensus set sizes, p-values, threshold comparisons, and gene
removal patterns across multiple runs. \# Create full convergence plot
\# conv_plot \<- full_gatai_convergence_plot(phyex_set, gatai_runs)
