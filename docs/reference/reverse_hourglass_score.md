# Reverse Hourglass Score Function

Compute the reverse hourglass score by comparing mid developmental
stages to early and late stages.

## Usage

``` r
reverse_hourglass_score(txi, modules)
```

## Arguments

- txi:

  Numeric vector of transcriptomic index values

- modules:

  A named list with elements 'early', 'mid', and 'late' containing stage
  indices for each developmental module

## Value

A numeric value representing the reverse hourglass score

## Details

The score is computed as the minimum of: - D1: mean(mid) - mean(early) -
D2: mean(mid) - mean(late)

Higher scores indicate stronger reverse hourglass patterns (mid stages
dominated by younger genes). \# Compute reverse hourglass score \#
modules \<- list(early = 1:3, mid = 4:6, late = 7:9) \# score \<-
reverse_hourglass_score(txi_values, modules)
