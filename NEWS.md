myTAI 0.0.2
===========

## New Features

- `combinatorialSignificance()`, `FlatLineTest()`, `ReductiveHourglassTest()`, and `EarlyConservationTest()`
now support multicore processing

- `MatchMap()` has been entirely rewritten and is now based on [dplyr](http://cran.rstudio.com/web/packages/dplyr/); additionally
it now has a new argument `accumulate` that allows you to accumulate multiple expression levels to a unique expressiion level for a unique gene id

## Vignettes

All three Vignettes: `Introduction`, `Intermediate`, and `Advanced` have been updated and extended.

## Bug Fixes

- two small bugs in `ReductiveHourglassTest()` and `EarlyConservationTest()` have been fixed that caused
that instead of displaying 3 or 4 plots (`par(mfrow=c(1,3))` or `par(mfrow=c(2,2))`) only 1 plot has been generated

- a small bug in `PlotMeans()` that caused the visualization of a wrong y-axis label
when plotting only one group of Phylostrata or Divergence Strata


myTAI 0.0.1
===========

Introducing myTAI 0.0.1:

A framework to perform phylotranscriptomics analyses for Evolutionary Developmental Biology research.