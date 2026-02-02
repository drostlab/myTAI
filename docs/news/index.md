# Changelog

## myTAI 2.3.5

- added functions `get_strata_legend` and `get_phylomap`
- removed faulty dynamic links from documentation

## myTAI 2.3.4

CRAN release: 2025-11-11

Patches: - taxid now throws warning instead of error and returns NULL if
no internet is available - taxonomy.Rmd vignetter now builds even with
no internet - fixed plot_contribution failing when some of the stratas
have no genes

## myTAI 2.3.0.9006

This update has focused on improving the single-cell phyloset
functionality.

The single-cell phylo expression object no longer depends on Seurat. You
can construct the ScPhyloExpressionSet either from a count matrix
(sparse or dense), using `ScPhyloExpressionSet_from_matrix`, or from a
Seurat object, using `ScPhyloExpressionSet_from_seurat`. For
consistency, one can use `BulkPhyloExpressionSet_from_df` instead of
`as_BulkPhyloExpressionSet`.

One key functionality of the single-cell object is the ability to switch
between different identities when plotting (equivalent to the
[`Seurat::Idents`](https://satijalab.github.io/seurat-object/reference/Idents.html)
functionality). This is done by setting the `::identities_label`
property of the object. The `::available_idents` property can be used to
see what options the user has in setting the current identity. By
setting `::idents_colours[[]]`, the user can choose a colour pallette
for the different identities when plotting, which are saved across
different plotting calls.

The computation of TAI values for single cell is now cached. Moreover,
we have readded the C++ accelerated code for the computation of TAI,
which upon profiling shows to be faster than the R version when handling
more than 100000 cells (an adaptive function chooses the appropriate
implementation for the object size).

Some of the plotting functionality had been improved and more options
were added for plotting (e.g. plot_gene_heatmap now allows for passing a
custom colour mapping for the rows (genes), instead of colouring them by
their strata; plot_signature for single cell should be more readable).

Bug fixes: - validation of S7 objects now works properly - printing of
object information now works properly (instead of dumping all the
properties information)

## myTAI 2.2.0.9006

Most importantly, new vignettes were added and the website has been
updated. Many thanks to
[@LotharukpongJS](https://github.com/LotharukpongJS)

Functionality wise:

- added more input validation across functions
- improved plotting in `destroy_pattern`
- renamed `tfPS` to `tf_PS` and fixed bug which prevented strata
  transformations from happening (PhyloExpressionSet now has an extra
  field [@strata_values](https://github.com/strata_values) which keeps
  track of the phylostratum values)
- added back standard deviation ribbon to plot_signature
- finer control over cell identity selection in the sc object
- fixed bugs in some of the plotting functions, such as
  plot_signature_multiple showing the colours in reverse,
  plot_genes_heatmap not working when given just one gene etc.

## myTAI 2.1.0.9000

### Major Refactoring and New Features

- **Single-cell support:** New `ScPhyloExpressionSet` S7 class for
  single-cell data, alongside `BulkPhyloExpressionSet` for bulk data.
- **Property renaming:** `conditions` → `identities` (and
  `conditions_label` → `identities_label`), `counts`/`counts_collapsed`
  → `expression`/`expression_collapsed`.
- **All core logic and plotting now robust to both bulk and single-cell
  objects.**

### Function Renaming and Prefixes

- **Statistical tests:** All core test functions now use the `stat_`
  prefix:
  - `flatline_test()` →
    [`stat_flatline_test()`](https://drostlab.github.io/myTAI/reference/stat_flatline_test.md)
  - `early_conservation_test()` →
    [`stat_early_conservation_test()`](https://drostlab.github.io/myTAI/reference/stat_early_conservation_test.md)
  - `late_conservation_test()` →
    [`stat_late_conservation_test()`](https://drostlab.github.io/myTAI/reference/stat_late_conservation_test.md)
  - `reductive_hourglass_test()` →
    [`stat_reductive_hourglass_test()`](https://drostlab.github.io/myTAI/reference/stat_reductive_hourglass_test.md)
  - `reverse_hourglass_test()` →
    [`stat_reverse_hourglass_test()`](https://drostlab.github.io/myTAI/reference/stat_reverse_hourglass_test.md)
  - `pairwise_test()` →
    [`stat_pairwise_test()`](https://drostlab.github.io/myTAI/reference/stat_pairwise_test.md)
  - `generic_conservation_test()` →
    [`stat_generic_conservation_test()`](https://drostlab.github.io/myTAI/reference/stat_generic_conservation_test.md)
  - `generate_conservation_txis()` →
    [`stat_generate_conservation_txis()`](https://drostlab.github.io/myTAI/reference/stat_generate_conservation_txis.md)
- **Gene selection/filtering:** All gene selection/filtering functions
  now use the `genes_` prefix:
  - `top_expression_genes()` →
    [`genes_top_mean()`](https://drostlab.github.io/myTAI/reference/genes_top_mean.md)
  - `top_variance_genes()` →
    [`genes_top_variance()`](https://drostlab.github.io/myTAI/reference/genes_top_variance.md)
  - `lowly_expressed_genes()` →
    [`genes_lowly_expressed()`](https://drostlab.github.io/myTAI/reference/genes_lowly_expressed.md)
  - `filter_dyn_expr()` →
    [`genes_filter_dynamic()`](https://drostlab.github.io/myTAI/reference/genes_filter_dynamic.md)

### Other Changes

- Many internal and Rd files renamed for consistency (e.g.,
  `gene_patterns.R` → `genes_patterns.R`, `S7_utils.R` → `utils_S7.R`).
- Deprecated/legacy files removed: `expression_utils.R`,
  `single_cell.R`.
- All documentation and examples updated to use new function/property
  names.

## myTAI 2.0.0.9000

### Major Changes

#### New S7 Class System

- Migrated from traditional R data structures to modern S7 classes
- New `PhyloExpressionSet` S7 class replaces the old data.frame-based
  format
- New `TestResult` S7 class for standardized storage of statistical test
  results
- Computed properties automatically calculate derived values like TXI,
  pTXI, conditions, etc.
- Built-in data validation and type checking through S7 properties
- Improved replicate data handling with automatic collapsing

#### Function Name Modernization

Function names have been updated to use snake_case convention:

##### Statistical Tests (old → new)

- `FlatLineTest()` → `flatline_test()`
- `ReductiveHourglassTest()` → `reductive_hourglass_test()`
- `EarlyConservationTest()` → `early_conservation_test()`
- `LateConservationTest()` → `late_conservation_test()`
- `ReverseHourglassTest()` → `reverse_hourglass_test()`
- `PairwiseTest()` → `pairwise_test()`

##### Visualization Functions (old → new)

- `PlotSignature()` →
  [`plot_signature()`](https://drostlab.github.io/myTAI/reference/plot_signature.md)
- `PlotPattern()` →
  [`plot_signature()`](https://drostlab.github.io/myTAI/reference/plot_signature.md)
- `PlotContribution()` →
  [`plot_contribution()`](https://drostlab.github.io/myTAI/reference/plot_contribution.md)
- `PlotDistribution()` →
  [`plot_distribution_strata()`](https://drostlab.github.io/myTAI/reference/plot_distribution_strata.md)
- `PlotCategoryExpr()` →
  [`plot_strata_expression()`](https://drostlab.github.io/myTAI/reference/plot_strata_expression.md)
- `PlotRE()` →
  [`plot_relative_expression_line()`](https://drostlab.github.io/myTAI/reference/plot_relative_expression_line.md)
- `PlotBarRE()` →
  [`plot_relative_expression_bar()`](https://drostlab.github.io/myTAI/reference/plot_relative_expression_bar.md)
- `PlotGeneSet()` →
  [`plot_gene_profiles()`](https://drostlab.github.io/myTAI/reference/plot_gene_profiles.md)
- `PlotMeans`, `PlotVars`, `PlotMedians` →
  \`plot_strata_expression(aggregate_FUN=“mean”/“var”/“median”)
- `PlotSignatureMultiple()` →
  [`plot_signature_multiple()`](https://drostlab.github.io/myTAI/reference/plot_signature_multiple.md)
- `PlotSignatureTransformed()` →
  [`plot_signature_transformed()`](https://drostlab.github.io/myTAI/reference/plot_signature_transformed.md)
- `PlotSignatureGeneQuantiles()` →
  [`plot_signature_gene_quantiles()`](https://drostlab.github.io/myTAI/reference/plot_signature_gene_quantiles.md)

##### Utility Functions (old → new)

- [`TAI()`](https://drostlab.github.io/myTAI/reference/TAI.md) →
  Computed property of PhyloExpressionSet (still accessible via
  [`TAI()`](https://drostlab.github.io/myTAI/reference/TAI.md))
- [`TDI()`](https://drostlab.github.io/myTAI/reference/TDI.md) →
  Computed property of PhyloExpressionSet (still accessible via
  [`TDI()`](https://drostlab.github.io/myTAI/reference/TDI.md))
- [`TEI()`](https://drostlab.github.io/myTAI/reference/TEI.md) →
  Computed property of PhyloExpressionSet (still accessible via
  [`TEI()`](https://drostlab.github.io/myTAI/reference/TEI.md))
- [`TPI()`](https://drostlab.github.io/myTAI/reference/TPI.md) →
  Computed property of PhyloExpressionSet (still accessible via
  [`TPI()`](https://drostlab.github.io/myTAI/reference/TPI.md))
- `pTAI()`, `pTDI` →
  [`sTXI()`](https://drostlab.github.io/myTAI/reference/sTXI.md)
  (generalized for all transcriptomic indices)
- `pMatrix()` →
  [`pTXI()`](https://drostlab.github.io/myTAI/reference/pTXI.md)
- `CollapseReplicates()` → `collapse`, Built into PhyloExpressionSet
  constructor
- `Expressed()` → `lowly_expressed_genes()`
- `MatchMap()` →
  [`match_map()`](https://drostlab.github.io/myTAI/reference/match_map.md)
- `SelectGeneSet()` →
  [`select_genes()`](https://drostlab.github.io/myTAI/reference/select_genes.md)
- `TopExpressionGenes()` → `top_expression_genes()`
- `TopVarianceGenes()` → `top_variance_genes()`
- `REMatrix()` →
  [`rel_exp_matrix()`](https://drostlab.github.io/myTAI/reference/rel_exp_matrix.md)
- `RE()` →
  [`relative_expression()`](https://drostlab.github.io/myTAI/reference/relative_expression.md)
- `omitMatrix()` →
  [`omit_matrix()`](https://drostlab.github.io/myTAI/reference/omit_matrix.md)
- `is.ExpressionSet()` → Built into S7 validation
- [`age.apply()`](https://drostlab.github.io/myTAI/reference/age.apply.md)
  →
  [`age.apply()`](https://drostlab.github.io/myTAI/reference/age.apply.md)
  (unchanged)
- [`tf()`](https://drostlab.github.io/myTAI/reference/tf.md) →
  [`tf()`](https://drostlab.github.io/myTAI/reference/tf.md) or
  [`transform_counts()`](https://drostlab.github.io/myTAI/reference/transform_counts.md)
- `tfPS()` → `tfPS()` (unchanged)
- `tfStability()` →
  [`tf_stability()`](https://drostlab.github.io/myTAI/reference/tf_stability.md)
- [`taxid()`](https://drostlab.github.io/myTAI/reference/taxid.md) →
  [`taxid()`](https://drostlab.github.io/myTAI/reference/taxid.md)
  (unchanged)

#### New Functions

- [`destroy_pattern()`](https://drostlab.github.io/myTAI/reference/destroy_pattern.md):
  Apply GATAI algorithm to identify pattern-contributing genes
- [`plot_signature_multiple()`](https://drostlab.github.io/myTAI/reference/plot_signature_multiple.md):
  Plot multiple signatures on the same plot
- [`plot_signature_gene_quantiles()`](https://drostlab.github.io/myTAI/reference/plot_signature_gene_quantiles.md):
  Plot signature with gene expression quantiles
- [`plot_signature_transformed()`](https://drostlab.github.io/myTAI/reference/plot_signature_transformed.md):
  Plot signatures with different transformations
- [`plot_sample_space()`](https://drostlab.github.io/myTAI/reference/plot_sample_space.md):
  Visualize sample relationships using PCA or UMAP
- [`plot_mean_var()`](https://drostlab.github.io/myTAI/reference/plot_mean_var.md):
  Mean-variance relationship plots
- [`plot_gene_profiles()`](https://drostlab.github.io/myTAI/reference/plot_gene_profiles.md):
  Individual gene expression profiles
- [`plot_distribution_expression()`](https://drostlab.github.io/myTAI/reference/plot_distribution_expression.md):
  Expression distribution plots
- [`plot_distribution_pTAI()`](https://drostlab.github.io/myTAI/reference/plot_distribution_pTAI.md):
  Partial TXI distribution plots
- [`plot_distribution_pTAI_qqplot()`](https://drostlab.github.io/myTAI/reference/plot_distribution_pTAI_qqplot.md):
  Q-Q plots for pTXI distributions
- [`plot_distribution_strata()`](https://drostlab.github.io/myTAI/reference/plot_distribution_strata.md):
  Phylostrata distribution plots
- [`plot_gene_heatmap()`](https://drostlab.github.io/myTAI/reference/plot_gene_heatmap.md):
  Gene expression heatmaps
- [`plot_gene_space()`](https://drostlab.github.io/myTAI/reference/plot_gene_space.md):
  Gene space visualization
- [`plot_cullen_frey()`](https://drostlab.github.io/myTAI/reference/plot_cullen_frey.md):
  Cullen-Frey plots for distribution assessment
- [`plot_null_txi_sample()`](https://drostlab.github.io/myTAI/reference/plot_null_txi_sample.md):
  Null TXI sample plots
- `as_PhyloExpressionSet()`: Convert data to PhyloExpressionSet S7
  object
- `get_sc_TAI()`: Single-cell TAI computation for Seurat objects=
- [`diagnose_test_robustness()`](https://drostlab.github.io/myTAI/reference/diagnose_test_robustness.md):
  Diagnose statistical test robustness
- [`remove_genes()`](https://drostlab.github.io/myTAI/reference/remove_genes.md):
  Remove genes from PhyloExpressionSet
- [`PS_colours()`](https://drostlab.github.io/myTAI/reference/PS_colours.md):
  Generate phylostratum color palettes
- [`ConservationTestResult()`](https://drostlab.github.io/myTAI/reference/ConservationTestResult.md):
  S7 class for conservation test results
- [`TestResult()`](https://drostlab.github.io/myTAI/reference/TestResult.md):
  S7 class for statistical test results

#### Enhanced Features

- Improved performance with parallelized C++ functions using
  RcppArmadillo
- All plots now use ggplot2 with consistent styling
- Comprehensive unit tests
- Updated to use modern R packages (ggplot2, dplyr, etc.)

#### Breaking Changes

- Function names have been updated to snake_case convention
- PhyloExpressionSet is now an S7 object instead of a data.frame
- Function signatures have been updated for consistency
- Some deprecated functions have been removed
- New package dependencies: S7, RcppArmadillo

#### Migration Guide

To migrate from myTAI 1.x to 2.0:

1.  **Convert data format**:

    ``` r
    # Old format (data.frame)
    old_phyex <- PhyloExpressionSetExample

    # New format (S7 object)
    new_phyex <- as_PhyloExpressionSet(old_phyex)
    ```

2.  **Update function calls**:

    ``` r
    # Old syntax
    PlotSignature(phyex_set, TestStatistic = "FlatLineTest")

    # New syntax
    plot_signature(phyex_set, conservation_test = flatline_test)
    ```

3.  **Access computed properties**:

    ``` r
    # Old syntax
    tai_values <- TAI(phyex_set)

    # New syntax
    tai_values <- phyex_set@TXI
    # or
    tai_values <- TAI(phyex_set)
    ```

## myTAI 1.0.2.9000

### New Functions

- New function `tfPS()` : Perform transformation of phylostratum values,
  analogous to `PS()` which transforms expression levels. Currently,
  `tfPS()` supports quantile rank transformation.
- New function `PairwiseTest()` : Statistically evaluate the pairwise
  difference in the phylotranscriptomic pattern between two contrasts
  based on or computations.
- New function `pairScore()` : Compute pairwise difference in TAI (or
  TDI) score.

### Bug and Issue Fixes

- Removed `taxonomy()` and references to it due to the deprecation of
  `taxize()`. This is needed for CRAN submission.

## myTAI 1.0.1.9000

### New Functions

- New function `PlotSignatureTransformed()` : Plot evolutionary
  signatures across transcriptomes and RNA-seq transformations

- New function `tfStability()`: Perform Permutation Tests Under
  Different Transformations (to test the robustness of the p-values of a
  given test (e.g. `FlatLineTest()`, `ReductiveHourglassTest()`,
  `ReverseHourglassTest()`, `EarlyConservationTest()` and
  `LateConservationTest()`) to expression data transformations)

- New function `LateConservationTest()` : Perform Reductive Late
  Conservation Test (to test for a high-mid-low (or high-high-low) TAI
  or TDI pattern)

- New function `lcScore()` : Compute the Hourglass Score for the
  LateConservationTest

- internal `rcpp` functions are now using Eigen to automatically enable
  parallelization

### New Features

- `PlotSignature()` updated to be able to perform the
  `TestStatistic = "LateConservationTest"`.

- `PlotSignature()` now prints p-value as a subtitle rather than via
  [`ggplot2::annotate()`](https://ggplot2.tidyverse.org/reference/annotate.html).

- [`tf()`](https://drostlab.github.io/myTAI/reference/tf.md) now has a
  `pseudocount` parameter, which is useful for performing logarithmic
  transformations when there are genes with 0 counts.

- [`tf()`](https://drostlab.github.io/myTAI/reference/tf.md) now
  supports `vst` and `rlog` transformations from
  [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

- [`tf()`](https://drostlab.github.io/myTAI/reference/tf.md) now has an
  `intergerise` parameter, which is needed when applying `vst` or `rlog`
  transformations.

- [`tf()`](https://drostlab.github.io/myTAI/reference/tf.md) updated
  documentation for performing rank transformation, which assigns ranks
  to the gene expression values within each stage, based on their
  relative positions compared to other values.

- Improvements to existing test functions (`ecScore()`, `rhScore()` and
  `reversehourglassScore()`) to give a message when the
  phylotranscriptomic pattern is unlikely to follow the test statistics.

- `FlatLineTest()` - newly returns the ks test statistics for the
  fitting of gamma

- `FlatLineTest()` - improved fitting

- `FlatLineTest()` - cpp functions are newly parallelized and progress
  bar is implemented for the computation of permutations

### Bug and Issue Fixes

- Some changes to remove errors and warnings from
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  and
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html),
  when building this package, which has been accumulated from previous
  updates.

## myTAI 1.0.0

### New Functions

- New function
  [`TEI()`](https://drostlab.github.io/myTAI/reference/TEI.md): Compute
  the Transcriptome Evolutionary Index
- New function `pMatrixTEI`: Compute Partial Transcriptome Evolutionary
  Index (TEI) Values
- New function `pStrataTEI`: Compute Partial Transcriptome Evolutionary
  Index (TEI) Strata Values
- New function `bootTEI`: Compute a Permutation Matrix of Transcriptome
  Evolutionary Index (TEI)
- new internal `rcpp` functions to support parallel C++ computations for
  [`TEI()`](https://drostlab.github.io/myTAI/reference/TEI.md),
  `pMatrixTEI()`, `pStrataTEI()`, `bootTEI()`

### New Features

- `CollapseReplicates()` now returns `tibble` objects
- `PlotCategoryExpr()` received a new argument `y.ticks`

### Bug and Issue Fixes

- `CollapseFromTo()` now has an exception when a replicate number `1` is
  passed to the function -\> previously this would cause an error to
  occur

## myTAI 0.9.3

CRAN release: 2021-02-24

- removing the depreciated `std::random_shuffle()` function to sample
  `plylostratum` or `divergence stratum` columns and replacing it with
  `std::shuffle()`. See full discussion
  [here](https://meetingcpp.com/blog/items/stdrandom_shuffle-is-deprecated.html).
- removing depreciated function calls such as
  [`dplyr::funs()`](https://dplyr.tidyverse.org/reference/funs.html) and
  [`tibble::is.tibble()`](https://tibble.tidyverse.org/reference/is.tibble.html)
- updated unit tests

## myTAI 0.9.2

CRAN release: 2020-01-10

- changing maintainer email-address

## myTAI 0.9.1

CRAN release: 2019-03-10

- fixing a unit test that uses `set.seed(123)` which causes an error in
  the new R version `3.6.0` due to the switch from a
  `non-uniform "Rounding" sampler` to a `"Rejection" sampler` in the new
  R version; the corresponding unit test `test-PlotEnrichment.R` was
  adjusted accordingly. Here the CRAN statement:

> Note that this ensures using the (old) non-uniform “Rounding” sampler
> for all 3.x versions of R, and does not add an R version dependency.
> Note also that the new “Rejection” sampler which R will use from 3.6.0
> onwards by default is definitely preferable over the old one, so that
> the above should really only be used as a temporary measure for
> reproduction of the previous behavior (and the run time tests relying
> on it).

## myTAI 0.9.0

CRAN release: 2019-02-06

#### New Functions

- new function `ReverseHourglassTest()` to perform a
  `Reverse Hourglass Test`. The Reverse Hourglass Test aims to
  statistically evaluate the existence of a reverse hourglass pattern
  based on TAI or TDI computations. The corresponding p-value quantifies
  the probability that a given TAI or TDI pattern (or any
  phylotranscriptomics pattern) does follow an hourglass like shape. A
  p-value \< 0.05 indicates that the corresponding phylotranscriptomics
  pattern does rather follow a reverse hourglass (low-high-low) shape.

- new function `reversehourglassScore()` for computing the
  `Reverse Hourglass Score` for the `Reverse Hourglass Test`

#### Updated Functionality

- function `PlotSignature()` receives a new `TestStatistic`
  (`TestStatistic = "ReverseHourglassTest"`) to perform a
  `revserse hourglass test` (= testing the significance of a
  low-high-low pattern)

## myTAI 0.8.0

CRAN release: 2018-05-23

#### Updated Functionality

- fix remaining issues when input is a `tibble`

## myTAI 0.7.0

CRAN release: 2018-04-11

### Updates

#### New Functions

- new function `PlotCIRatio()` to compute and visualize TAI/TDI etc
  patters using bootstrapping and confidence intervals (contributed by
  [@ljljolinq1010](https://github.com/ljljolinq1010))

#### Update Functionality

- all functions can now handle `tibble` data as input -\> before there
  were errors thrown when input data wasn’t in strict `data.frame`
  format

## myTAI 0.6.0

CRAN release: 2017-07-03

#### Updates

- `is.ExpressionSet()` now prints out more detailed error messages when
  ExpressionSet is violated

- adapt `PlotContribution()` to new version of `dplyr` where
  [`summarise_each()`](https://dplyr.tidyverse.org/reference/summarise_each.html)
  is deprecated.

Error message occurring after new `dplyr` release was:

2.  Failure: `PlotContribution()` works properly with
    DivergenceExpressionSet input…
    ([@test-PlotContribution](https://github.com/test-PlotContribution).R#16)
    PlotContribution(DivergenceExpressionSetExample, legendName = “DS”)
    produced messages.

[`summarise_each()`](https://dplyr.tidyverse.org/reference/summarise_each.html)
is deprecated. Use
[`summarise_all()`](https://dplyr.tidyverse.org/reference/summarise_all.html),
[`summarise_at()`](https://dplyr.tidyverse.org/reference/summarise_all.html)
or
[`summarise_if()`](https://dplyr.tidyverse.org/reference/summarise_all.html)
instead. To map `funs` over all variables, use
[`summarise_all()`](https://dplyr.tidyverse.org/reference/summarise_all.html)
[`summarise_each()`](https://dplyr.tidyverse.org/reference/summarise_each.html)
is deprecated. Use
[`summarise_all()`](https://dplyr.tidyverse.org/reference/summarise_all.html),
[`summarise_at()`](https://dplyr.tidyverse.org/reference/summarise_all.html)
or
[`summarise_if()`](https://dplyr.tidyverse.org/reference/summarise_all.html)
instead. To map `funs` over all variables, use
[`summarise_all()`](https://dplyr.tidyverse.org/reference/summarise_all.html)

Is now fixed.

## myTAI 0.5.0

CRAN release: 2017-03-14

#### New Functions

- new function `PlotSignature()` allows users to plot evolutionary
  signatures across transcriptomes (based on ggplot2 -\> new main
  visualization function aiming to replace the `PlotPattern()` function)

- new function
  [`TPI()`](https://drostlab.github.io/myTAI/reference/TPI.md) allows
  users to compute the Transcriptome Polymorphism Index introduced by
  `Gossmann et al., 2015`.

- new function `PlotMedians()` allows users to compute and visualize the
  median expression of all age categories

- new function `PlotVars()` allows users to compute and visualize the
  expression variance of all age categories

#### Updates

- `PlotContribution()` is now based on ggplot2 and loses base graphics
  arguments

- now R/RcppExports.R and src/rcpp_funcs.cpp are included in the package
  due to previous compilation problems (see also [stackoverflow
  discussion](https://stackoverflow.com/questions/34585560/travis-ci-r-package-error-in-documentation))

- `MatchMap()` is now based on
  [`dplyr::inner_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
  to match age category table with a gene expression dataset

- `PlotCorrelation()` has been extended and optimized for producing high
  publication quality plots

- `PlotMeans()` is now based on ggplot2 and lost all base graphics
  arguments.

- `PlotRE()` is now based on ggplot2 and lost all base graphics
  arguments.

#### Vignettes

- In `Introduction` vignette: complete restructuring of the Introduction
- In `Introduction` vignette: add new ggplot2 based examples

## myTAI 0.4.0

CRAN release: 2016-05-04

#### New Functions

- a new function `PlotSelectedAgeDistr()` allowing unsers to visualize
  the PS or DS gene distribution of a subset of genes stored in the
  input ExpressionSet object
- a new function `PlotGroupDiffs()` allowing users to plot the
  significant differences between gene expression distributions of PS or
  DS groups
- a new function `GroupDiffs()` allowing users to perform statistical
  tests to quantify the gene expression level differences between all
  genes of defined PS or DS groups

#### Updates

- `PlotDistribution()` now uses ggplot2 to visualize the PS or DS
  distribution and is also based on the new function
  `PlotSelectedAgeDistr()`; furthermore it loses arguments `plotText`
  and `...` and gains a new argument `legendName`

- remove arguments ‘main.text’ and ‘…’ from `PlotCorrelation()`

- `PlotCorrelation()` is now based on ggplot2

- `PlotGroupDiffs()` receives a new argument `gene.set` allowing users
  to statistically quantify the group specific PS/DS differences of a
  selected set of genes

- analogously to `PlotGroupDiffs()` the function `GroupDiffs()` also
  receives a new argument `gene.set` allowing users to statistically
  quantify the group specific PS/DS differences of a selected set of
  genes

- Fixing wrong x-axis labeling in `PlotCategoryExpr()` when
  `type = "stage-centered"` is specified

- `PlotCategoryExpr()` now also prints out the PS/DS absolute frequency
  distribution of the selected `gene.set`

## myTAI 0.3.0

CRAN release: 2015-07-31

#### Vignettes

- adding examples for `PlotCategoryExpr()` to `Advanced` Vignette
- adding examples for `PlotReplicateQuality()` to `Expression` vignette

#### New Functions

- a new function `PlotCategoryExpr()` allowing users to plot the
  expression levels of each age or divergence category as boxplot, dot
  plot or violin plot
- a new function `PlotReplicateQuality()` allowing users to visualize
  the quality of biological replicates

## myTAI 0.2.1

CRAN release: 2015-07-23

#### Vignettes

- fixed a wrong example in the Enrichment vignette
  (<https://github.com/HajkD/myTAI/commit/8d52fd60c274361dc9028dec3409abf60a738d8a>)

#### Updates

- `PlotGeneSet()` and `SelectGeneSet()` now have a new argument
  `use.only.map` specifying whether or not instead of using a standard
  `ExpressionSet` a `Phylostratigraphic Map` or `Divergene Map` is
  passed to the function.
- a wrong version of the edgeR Bioconductor package was imported causing
  version 0.2.0 to fail R CMD Check on unix based systems

## myTAI 0.2.0

CRAN release: 2015-07-10

#### Vignettes

- adding new vignette **Taxonomy** providing spep by step instructions
  on retrieving taxonomic information for any organism of interest

- adding new vignette **Expression Analysis** providing use cases to
  perform gene expression data analysis with `myTAI`

- adding new vignette **Enrichment** providing step-by-step instructions
  on how to perform PS and DS enrichment analyses with
  `PlotEnrichment()`

- adding examples for `pStrata()`, `pMatrix()`, `pTAI()`, `pTDI()`, and
  `PlotContribution()` to the **Introduction** Vignette

#### New Functions

- a new function `taxonomy()` allows users to retrieve taxonomic
  information for any organism of interest; this function has been taken
  from the [biomartr](https://github.com/ropensci/biomartr) package and
  was removed from `biomartr` afterwards. Please notice, that in myTAI
  version 0.1.0 the Introduction vignette referenced to the `taxonomy()`
  function in `biomartr`. This is no longer the case (since myTAI
  version 0.2.0), because now `taxonomy()` is implemented in myTAI.

- the new `taxonomy()` function is based on the powerful R package
  [taxize](https://github.com/ropensci/taxize).

- a new function `SelectGeneSet()` allows users to fastly select a
  subset of genes in an ExpressionSet

- a new function `DiffGenes()` allows users to perform differential gene
  expression analysis with ExpressionSet objects

- a new function `EnrichmentTest()` allows users to perform a Fisher’s
  exact test based enrichment analysis of over or underrepresented
  Phylostrata or Divergence Strata within a given gene set without
  having to plot the result

- a new function `PlotGeneSet()` allows users to visualize the
  expression profiles of a given gene set

- a new function `PlotEnrichment()` allows users to visualize the
  Phylostratum or Divergence Stratum enrichment of a given Gene Set as
  well as computing Fisher’s exact test to quantify the statistical
  significance of enrichment

- a new function `PlotContribution()` allows users to visualize the
  Phylostratum or Divergence Stratum contribution to the global TAI/TDI
  pattern

- a new function `pTAI()` allows users to compute the phylostratum
  contribution to the global TAI pattern

- a new function `pTDI()` allows users to compute the divergence stratum
  contribution to the global TDI pattern

#### Updates

- `FilterRNASeqCT()` has been renamed to `Expressed()` allowing users to
  apply this filter function to RNA-Seq data as well as to microarray
  data
- `PlotRE()` and `PlotMeans()` are now based on colors from the
  RColorBrewer package (default)
- `PlotRE()` and `PlotMeans()` now have a new argument `colors` allowing
  unsers to choose custom colors for the visualized relative or mean
  expression profiles
- `geom.mean()` and `harm.mean()` now are external functions accessible
  to the `myTAI` user

## myTAI 0.1.0

CRAN release: 2015-05-24

#### Main News

- now all functions have unit tests

#### New Functions

- a new function `pStrata()` allows users to compute partial TAI/TDI
  values for all Phylostrata or Divergence Strata

- a new function `CollapseReplicates()` allows users to combine
  replicate expression levels in ExpressionSet objects

- a new function `FilterRNASeqCT()` allows users to filter expression
  levels of `ExpressionSet` objects deriving from RNA-Seq count tables

#### Updates

- function `MatchMap()` now receives a new argument `remove.duplicates`
  allowing users to delete duplicate gene ids (that might be stored in
  the input PhyoMap or DivergenceMap) during the process of matching a
  Map with an ExpressionSet

- `FlatLineTest()`, `ReductiveHourglassTest()`,
  `EarlyConservationTest()`, and `PlotPattern()` implement a new
  argument `custom.perm.matrix` allowing users to pass their own
  (custom) permutation matrix to the corresponding function. All
  subsequent test statistics and p-value/std.dev computations are then
  based on this custom permutation matrix

- `EarlyConservationTest()` and `ReductiveHourglassTest()` now have a
  new parameter `gof.warning` allowing users to choose whether or not
  non significant goodness of fit results should be printed as warning

- now when specifying `TestStatistic = NULL` in `PlotPattern()` only the
  TAI/TDI profile is drawn (without performing any test statistics);
  this is equavalent to performing:
  `plot(TAI(PhyloExpressionSetExample)`

- function `combinatorialSignificance()` is now named
  `CombinatorialSignificance()`

- changing the title and description of the `myTAI` package

- some minor changes in vignettes and within the documentation of
  functions

## myTAI 0.0.2

CRAN release: 2014-12-27

### New Features in v. 0.0.2

- `combinatorialSignificance()`, `FlatLineTest()`,
  `ReductiveHourglassTest()`, and `EarlyConservationTest()` now support
  multicore processing

- `MatchMap()` has been entirely rewritten and is now based on
  [dplyr](https://github.com/tidyverse/dplyr); additionally it now has a
  new argument `accumulate` that allows you to accumulate multiple
  expression levels to a unique expressiion level for a unique gene id

#### Vignettes

All three Vignettes: `Introduction`, `Intermediate`, and `Advanced` have
been updated and extended.

#### Bug Fixes

- two small bugs in `ReductiveHourglassTest()` and
  `EarlyConservationTest()` have been fixed that caused that instead of
  displaying 3 or 4 plots (`par(mfrow=c(1,3))` or `par(mfrow=c(2,2))`)
  only 1 plot has been generated

- a small bug in `PlotMeans()` that caused the visualization of a wrong
  y-axis label when plotting only one group of Phylostrata or Divergence
  Strata

## myTAI 0.0.1

CRAN release: 2014-11-01

Introducing myTAI 0.0.1:

A framework to perform phylotranscriptomics analyses for Evolutionary
Developmental Biology research.
