# Create Convergence Plots for GATAI Analysis

Create plots showing how consensus gene set sizes and p-values converge
across GATAI runs for different threshold values.

## Usage

``` r
convergence_plots(phyex_set, runs, ps = c(0.5))
```

## Arguments

- phyex_set:

  A PhyloExpressionSet object

- runs:

  List of GATAI run results

- ps:

  Vector of consensus thresholds to analyze (default: c(0.5))

## Value

A list with two ggplot objects:

- counts:

  Plot showing convergence of consensus set sizes

- pval:

  Plot showing convergence of p-values

## Details

This function analyzes how consensus gene sets and their statistical
significance change as more GATAI runs are included in the analysis. It
uses cached null distributions for efficient p-value calculation. \#
Create convergence plots \# conv_plots \<- convergence_plots(phyex_set,
gatai_runs, ps = c(0.25, 0.5, 0.75))
