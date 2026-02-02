# Early Conservation Score Function

Compute the early conservation score by comparing mid and late
developmental stages to early stages.

## Usage

``` r
ec_score(txi, modules)
```

## Arguments

- txi:

  Numeric vector of transcriptomic index values

- modules:

  A named list with elements 'early', 'mid', and 'late' containing stage
  indices for each developmental module

## Value

A numeric value representing the early conservation score

## Details

The score is computed as the minimum of: - D1: mean(mid) - mean(early) -
D2: mean(late) - mean(early)

Higher scores indicate stronger early conservation patterns. \# Compute
early conservation score \# modules \<- list(early = 1:3, mid = 4:6,
late = 7:9) \# score \<- ec_score(txi_values, modules)
