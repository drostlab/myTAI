# Pairwise Score Function

Compute the pairwise contrast score between two groups of developmental
stages.

## Usage

``` r
pair_score(txi, modules, alternative = c("greater", "less"))
```

## Arguments

- txi:

  Numeric vector of transcriptomic index values

- modules:

  A named list with elements 'contrast1' and 'contrast2' containing
  stage indices for each contrast group

- alternative:

  Character string specifying the alternative hypothesis: "greater" or
  "less"

## Value

A numeric value representing the pairwise contrast score

## Details

The score is computed as mean(contrast1) - mean(contrast2). For
alternative = "less", the score is negated. \# Compute pairwise score \#
modules \<- list(contrast1 = 1:3, contrast2 = 7:9) \# score \<-
pair_score(txi_values, modules, "greater")
