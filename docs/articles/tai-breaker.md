# Break TAI patterns using gaTAI

## Motivation for breaking the hourglass

Inspired by the ‘ontogeny-phylogeny’ correlations and the hourglass
model of evolutionary development, we recently asked whether we could
break the hourglass or any TAI profile by extracting a stable set of
genes, and gain information about the role these extracted genes have in
the development of the organism.

Further motivation is in [this
thesis](https://dspace.cuni.cz/handle/20.500.11956/193428).

## Break the hourglass using gaTAI

``` r
library(myTAI)
data("example_phyex_set")
```

``` r
myTAI::destroy_pattern(example_phyex_set)
```

…

This is a teaser.

The
[`myTAI::destroy_pattern()`](https://drostlab.github.io/myTAI/reference/destroy_pattern.md)
requires the package `gataiR` which is still under development.
