# Reductive Hourglass Score Function

Compute the reductive hourglass score by comparing early and late
developmental stages to mid stages.

## Usage

``` r
reductive_hourglass_score(txi, modules)
```

## Arguments

- txi:

  Numeric vector of transcriptomic index values

- modules:

  A named list with elements 'early', 'mid', and 'late' containing stage
  indices for each developmental module

## Value

A numeric value representing the reductive hourglass score

## Details

The score is computed as the minimum of: - D1: mean(early) - mean(mid) -
D2: mean(late) - mean(mid)

Higher scores indicate stronger reductive hourglass patterns (mid stages
dominated by older genes). \# Compute reductive hourglass score \#
modules \<- list(early = 1:3, mid = 4:6, late = 7:9) \# score \<-
reductive_hourglass_score(txi_values, modules)
