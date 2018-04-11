myTAI 0.7.0
===========

## Updates

### New Functions

- new function `PlotCIRatio()` to compute and visualize TAI/TDI etc patters using bootstrapping and confidence intervals (contributed by @ljljolinq1010)


### Update Functionality

- all functions can now handle `tibble` data as input -> before there were errors
thrown when input data wasn't in strict `data.frame` format

myTAI 0.6.0
===========

### Updates

- `is.ExpressionSet()` now prints out more detailed error messages when [ExpressionSet format](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd#defining-input-data-standards) is violated

- adapt `PlotContribution()` to new version of `dplyr` where `summarise_each()` is deprecated.

Error message accuring after new `dplyr` release was:

2. Failure: `PlotContribution()` works properly with DivergenceExpressionSet input... (@test-PlotContribution.R#16) 
  PlotContribution(DivergenceExpressionSetExample, legendName = "DS") produced messages.
  
  
  `summarise_each()` is deprecated.
  Use `summarise_all()`, `summarise_at()` or `summarise_if()` instead.
  To map `funs` over all variables, use `summarise_all()`
  `summarise_each()` is deprecated.
  Use `summarise_all()`, `summarise_at()` or `summarise_if()` instead.
  To map `funs` over all variables, use `summarise_all()`
  
Is now fixed.  


myTAI 0.5.0
===========

### New Functions

- new function `PlotSignature()` allows users to plot evolutionary signatures across transcriptomes (based on ggplot2 -> new main visualization function aiming to replace the `PlotPattern()` function)

- new function `TPI()` allows users to compute the Transcriptome Polymorphism Index introduced by `Gossmann et al., 2015`.

- new function `PlotMedians()` allows users to compute and visualize the median expression of all age categories

- new function `PlotVars()` allows users to compute and visualize the expression variance of all age categories

### Updates

- `PlotContribution()` is now based on ggplot2 and loses base graphics arguments

- now R/RcppExports.R and src/rcpp_funcs.cpp are included in the package due to previous compilation problems (see also [stackoverflow discussion](http://stackoverflow.com/questions/34585560/travis-ci-r-package-error-in-documentation))

- `MatchMap()` is now based on `dplyr::inner_join()` to match age category table with a gene expression dataset

- `PlotCorrelation()` has been extended and optimized for producing high publication quality plots 

- `PlotMeans()` is now based on ggplot2 and lost all base graphics arguments.

- `PlotRE()` is now based on ggplot2 and lost all base graphics arguments.

### Vignettes

- In `Introduction` vignette: complete restructuring of the Introduction 
- In `Introduction` vignette: add new ggplot2 based examples

 
 
myTAI 0.4.0
===========

### New Functions

- a new function `PlotSelectedAgeDistr()` allowing unsers to visualize the PS or DS gene distribution of a subset of genes stored in the input ExpressionSet object
- a new function `PlotGroupDiffs()` allowing users to plot the significant differences between gene expression distributions of PS or DS groups
- a new function `GroupDiffs()` allowing users to perform statistical tests to quantify the gene expression level differences between all genes of defined PS or DS groups 

### Updates

- `PlotDistribution()` now uses ggplot2 to visualize the PS or DS distribution and
is also based on the new function `PlotSelectedAgeDistr()`; furthermore it loses arguments `plotText` and `...` and gains a new argument `legendName`

- remove arguments 'main.text' and '...' from `PlotCorrelation()`
- `PlotCorrelation()` is now based on ggplot2
- `PlotGroupDiffs()` receives a new argument `gene.set` allowing users to statistically quantify the group specific PS/DS differences of a selected set of genes
- analogously to `PlotGroupDiffs()` the function `GroupDiffs()` also receives a new argument `gene.set` allowing users to statistically quantify the group specific PS/DS differences of a selected set of genes
- Fixing wrong x-axis labeling in `PlotCategoryExpr()` when `type = "stage-centered"` is specified
- `PlotCategoryExpr()` now also prints out the PS/DS absolute frequency distribution of the selected `gene.set`

myTAI 0.3.0
===========

### Vignettes

- adding examples for `PlotCategoryExpr()` to `Advanced` Vignette
- adding examples for `PlotReplicateQuality()` to `Expression`
 vignette
 
### New Functions

- a new function `PlotCategoryExpr()` allowing users to plot the expression levels of each age or divergence category as boxplot, dot plot or violin plot
- a new function `PlotReplicateQuality()` allowing users to visualize the quality of biological replicates


myTAI 0.2.1
===========

### Vignettes

- fixed a wrong example in the Enrichment vignette (https://github.com/HajkD/myTAI/commit/8d52fd60c274361dc9028dec3409abf60a738d8a)


### Updates

- `PlotGeneSet()` and `SelectGeneSet()` now have a new argument `use.only.map` specifying whether or not instead of using a standard `ExpressionSet` a `Phylostratigraphic Map` or `Divergene Map` is passed to the function.
- a wrong version of the edgeR Bioconductor package was imported causing version 0.2.0 to fail R CMD Check on unix based systems


myTAI 0.2.0
===========

### Vignettes

- adding new vignette __Taxonomy__ providing spep by step instructions on retrieving taxonomic information for any organism of interest

- adding new vignette __Expression Analysis__ providing use cases to perform gene expression
data analysis with `myTAI`

- adding new vignette __Enrichment__ providing step-by-step instructions on how to perform PS and DS enrichment analyses with `PlotEnrichment()`

- adding examples for `pStrata()`, `pMatrix()`, `pTAI()`, `pTDI()`, and `PlotContribution()`
to the __Introduction__ Vignette

### New Functions

- a new function `taxonomy()` allows users to retrieve taxonomic information for any organism of interest; this function has been taken from the [biomartr](https://github.com/HajkD/biomartr) package and was removed from `biomartr` afterwards. Please notice, that in myTAI version 0.1.0 the Introduction vignette referenced to the `taxonomy()` function in `biomartr`. This is no longer the case (since myTAI version 0.2.0), because now `taxonomy()` is implemented in myTAI. 

- the new `taxonomy()` function is based on the powerful R package [taxize](https://github.com/ropensci/taxize).

- a new function `SelectGeneSet()` allows users to fastly select a subset of genes in an ExpressionSet

- a new function `DiffGenes()` allows users to perform differential gene expression analysis with ExpressionSet objects

- a new function `EnrichmentTest()` allows users to perform a Fisher's exact test based enrichment analysis of over or underrepresented Phylostrata or Divergence Strata within a given gene set without having to plot the result

- a new function `PlotGeneSet()` allows users to visualize the expression profiles of a given gene set

- a new function `PlotEnrichment()` allows users to visualize the Phylostratum or Divergence Stratum enrichment of a given Gene Set as well as computing Fisher's exact test to quantify the statistical significance of enrichment

- a new function `PlotContribution()` allows users to visualize the Phylostratum or Divergence Stratum contribution to the global TAI/TDI pattern

- a new function `pTAI()` allows users to compute the phylostratum contribution to the global TAI pattern

- a new function `pTDI()` allows users to compute the divergence stratum contribution to the global TDI pattern

### Updates

- `FilterRNASeqCT()` has been renamed to `Expressed()` allowing users to apply this filter function
to RNA-Seq data as well as to microarray data 
- `PlotRE()` and `PlotMeans()` are now based on colors from the RColorBrewer package (default)
- `PlotRE()` and `PlotMeans()` now have a new argument `colors` allowing unsers to choose custom colors for the visualized relative or mean expression profiles 
- `geom.mean()` and `harm.mean()` now are external functions accessible to the `myTAI` user

myTAI 0.1.0
===========

### Main News

- now all functions have unit tests

### New Functions

- a new function `pStrata()` allows users to compute partial TAI/TDI values for all Phylostrata or Divergence Strata

- a new function `CollapseReplicates()` allows users to combine replicate expression levels in ExpressionSet objects

- a new function `FilterRNASeqCT()` allows users to filter expression levels of `ExpressionSet` objects deriving from RNA-Seq count tables  

### Updates

- function `MatchMap()` now receives a new argument `remove.duplicates` allowing users to delete
duplicate gene ids (that might be stored in the input PhyoMap or DivergenceMap) during the process of matching a Map with an ExpressionSet

- `FlatLineTest()`, `ReductiveHourglassTest()`, `EarlyConservationTest()`, and `PlotPattern()` implement a new argument `custom.perm.matrix` allowing users to pass their own (custom) permutation matrix to the corresponding function. All subsequent test statistics and p-value/std.dev computations are then based on this custom permutation matrix

- `EarlyConservationTest()` and `ReductiveHourglassTest()` now have a new parameter `gof.warning` allowing users to choose whether or not non significant goodness of fit results should be printed as warning

- now when specifying `TestStatistic = NULL` in `PlotPattern()` only the TAI/TDI profile is drawn (without performing any test statistics); this is equavalent to performing: `plot(TAI(PhyloExpressionSetExample)`

- function `combinatorialSignificance()` is now named `CombinatorialSignificance()`

- changing the title and description of the `myTAI` package 

- some minor changes in vignettes and within the documentation of functions


myTAI 0.0.2
===========

## New Features in v. 0.0.2

- `combinatorialSignificance()`, `FlatLineTest()`, `ReductiveHourglassTest()`, and `EarlyConservationTest()`
now support multicore processing

- `MatchMap()` has been entirely rewritten and is now based on [dplyr](https://github.com/hadley/dplyr); additionally
it now has a new argument `accumulate` that allows you to accumulate multiple expression levels to a unique expressiion level for a unique gene id

### Vignettes

All three Vignettes: `Introduction`, `Intermediate`, and `Advanced` have been updated and extended.

### Bug Fixes

- two small bugs in `ReductiveHourglassTest()` and `EarlyConservationTest()` have been fixed that caused
that instead of displaying 3 or 4 plots (`par(mfrow=c(1,3))` or `par(mfrow=c(2,2))`) only 1 plot has been generated

- a small bug in `PlotMeans()` that caused the visualization of a wrong y-axis label
when plotting only one group of Phylostrata or Divergence Strata


myTAI 0.0.1
===========

Introducing myTAI 0.0.1:

A framework to perform phylotranscriptomics analyses for Evolutionary Developmental Biology research.
