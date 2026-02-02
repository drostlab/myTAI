# Create Threshold Comparison Plots

Create plots comparing how different consensus thresholds affect gene
set sizes and p-values in GATAI analysis.

## Usage

``` r
threshold_comparison_plots(phyex_set, runs)
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- runs:

  List of GATAI run results

## Value

A list with two ggplot objects:

- counts:

  Plot showing gene counts across thresholds

- pval:

  Plot showing p-values across thresholds

## Details

This function analyzes how the choice of consensus threshold (minimum
number of runs a gene must appear in) affects the final gene set size
and statistical significance. Uses cached null distributions for
efficient p-value calculation. \# Create threshold comparison plots \#
thresh_plots \<- threshold_comparison_plots(phyex_set, gatai_runs)
