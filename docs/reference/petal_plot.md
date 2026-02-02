# Create Petal Plot for Gene Removal Analysis

Create a petal plot showing how many genes are removed per run relative
to the consensus set.

## Usage

``` r
petal_plot(sets, p = 0.5)
```

## Arguments

- sets:

  List of gene sets from GATAI runs

- p:

  Consensus threshold (default: 0.5)

## Value

A ggplot2 petal plot

## Details

This function creates a petal plot visualization showing the
relationship between individual GATAI runs and the consensus gene set,
highlighting how many genes are unique to each run. \# Create petal plot
\# petal \<- petal_plot(gatai_runs, p = 0.5)
