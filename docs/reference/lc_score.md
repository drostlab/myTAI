# Late Conservation Score Function

Compute the late conservation score by comparing early and mid
developmental stages to late stages.

## Usage

``` r
lc_score(txi, modules)
```

## Arguments

- txi:

  Numeric vector of transcriptomic index values

- modules:

  A named list with elements 'early', 'mid', and 'late' containing stage
  indices for each developmental module

## Value

A numeric value representing the late conservation score

## Details

The score is computed as the minimum of: - D1: mean(early) -
mean(late) - D2: mean(mid) - mean(late)

Higher scores indicate stronger late conservation patterns. \# Compute
late conservation score \# modules \<- list(early = 1:3, mid = 4:6, late
= 7:9) \# score \<- lc_score(txi_values, modules)
