# Other evolutionary and expression indices

## Motivation for non-TAI approaches

So far, most of the vignettes have focused on workflows and discussions
centred around the transcriptome age index (TAI). However, we can also
assign other information for each gene (analogous to gene age in TAI).
Examples include [`TDI`](https://www.nature.com/articles/nature11394),
[`TSI`](https://www.nature.com/articles/s41586-024-08059-8),
[`TPI`](https://doi.org/10.1093/molbev/msw044) and so on (though the
names of these indices can differ between studies and are thus not
standardised).

In the sections below, we describe the workflow for the transcriptome
divergence index (TDI) and transcriptome specificity index (TSI). These
and other index-based methods use the same framework as TAI, aside from
the method to give a factor or numeric value to each gene (e.g. gene age
for TAI).

The good news is that `myTAIv2` is rather flexible; both
`BulkPhyloExpressionSet` and `ScPhyloExpressionSet` objects can
accommodate other indices, along with most of the `myTAIv2` functions.

Users can simply follow the workflow outlined in [statistical
tests](https://drostlab.github.io/myTAI/articles/tai-stats.md),
[transformations](https://drostlab.github.io/myTAI/articles/tai-transform.md)
along with most of the [plotting
functionalities](https://drostlab.github.io/myTAI/articles/tai-gallery.md),
using non-TAI data.

## Transcriptome divergence index (TDI)

The `TDI` ([Quint et al.,
2012](https://www.nature.com/articles/nature11394)) uses deciled dNdS
values for each gene and summarises this information across all genes at
a given stage (analogous to TAI) to quantify the overall transcriptome
‘divergence’.

For the vast majority of genes, the dNdS values are below 1
(i.e. purifying selection). Therefore, while TDI is called transcriptome
*divergence* index, TDI can be thought to capture the degree of
purifying selection ([Lotharukpong et
al. 2024](https://www.nature.com/articles/s41586-024-08059-8)).

R packages such as [`orthologr`](https://github.com/drostlab/orthologr)
can facilitate dNdS estimation and
[divergence-stratigraphy](https://drostlab.github.io/orthologr/articles/divergence_stratigraphy.html)
computation used in the TDI approach.

Example TDI workflow using mock data

Assuming that we have the [divergence
strata](https://drostlab.github.io/orthologr/articles/divergence_stratigraphy.html),
similar to `TAI`, `TDI` is calculated as:

``` math
{TDI}_s = \frac{\sum_{i=1}^{n} {ds}_i \cdot e_{is}}{\sum_{i=1}^{n} e_{is}}
```

where $`e_{is}`$ is the expression level of gene $`i`$ at a given sample
$`s`$ (e.g. a biological replicate for a developmental stage), and
$`{ds}_i`$ is the evolutionary age of gene $`i`$.

Like TAI, for developmental time-course data, $`{TDI}_s`$ can be grouped
for each stage (if replicate data exists) and compared across
developmental stages to obtain the overall TDI profile.  
For pairwise comparisons (see
[📈](https://drostlab.github.io/myTAI/articles/tai-stats.md)),
$`{TDI}_s`$ can be grouped as one of the conditions to be compared.  
For single cell data, (see
[📊](https://drostlab.github.io/myTAI/articles/phylo-expression-object.md)),
$`{TDI}_s`$ is the TDI for a given cell or cell type (if pseudo-bulked).

### Mock bulk RNA-seq dataset for TDI

``` r
library(myTAI)
data("example_phyex_set_old")
```

Now we create the mock dataset, where we assign a value between 1 and 10
for all genes. This is the same distribution of values as the
`divergence strata`.

``` r
# this is a mock divergence strata dataset
# we assign a value between 1 and 10 for all genes
gene_names <- example_phyex_set_old@gene_ids
example_divergence_strata <- 
    setNames(sample(1:10, length(gene_names), replace = TRUE), gene_names) |>
    as.factor()
# mock expression data taken from example_phyex_set_old@expression
example_expression <- 
    example_phyex_set_old@expression |> 
    as.data.frame() |> 
    tibble::rownames_to_column(var = "GeneID")
```

Now we construct the input data.frame for
[`as_BulkPhyloExpressionSet()`](https://drostlab.github.io/myTAI/reference/as_BulkPhyloExpressionSet.md)

``` r
example_phyex_set_tdi.df <-
    data.frame(phylorank = example_divergence_strata) |>
    tibble::rownames_to_column(var = "GeneID") |>
    dplyr::select(phylorank, GeneID) |>
    dplyr::left_join(example_expression, by = "GeneID")
```

Create the `BulkPhyloExpressionSet`.

``` r
example_phyex_set_tdi <- 
    myTAI::as_BulkPhyloExpressionSet(
        example_phyex_set_tdi.df,
        name = "Example TDI data",
        index_type = "TDI")
```

plot away!

``` r
myTAI::plot_signature(example_phyex_set_tdi)
```

    ## Computing: [========================================] 100% (done)

![plot signature function output for mock TDI
data](other-strata_files/figure-html/unnamed-chunk-5-1.png)

Looks interesting right?

… although as you can expect from the
[`myTAI::stat_flatline_test()`](https://drostlab.github.io/myTAI/reference/stat_flatline_test.md),
because we assigned random values, the overall pattern doesn’t deviate
from a flat line.

``` r
myTAI::stat_flatline_test(example_phyex_set_tdi)
```

![stat_flatline_test function
output](other-strata_files/figure-html/unnamed-chunk-6-1.png)

The dNdS values are not without issues. For example, it is important to
choose an appropriate evolutionary distance to pairwise align and
compare nucleotide sequences. Too close results in no meaningful
substitution and too distant results in saturation. To overcome this,
similar measures to the dNdS can be computed from MSAs.

## Transcriptome specificity index (TSI)

The `TSI` ([Lotharukpong et al.,
2024](https://www.nature.com/articles/s41586-024-08059-8)) uses
expression breadth/specificity information that can be obtained for each
gene and summarises this information across all genes at a given stage
(analogous to TAI) to quantify the overall transcriptome ‘specificity’.

Similar to `TAI`, `TSI` is calculated as:

``` math
{TSI}_s = \frac{\sum_{i=1}^{n} {ts}_i \cdot e_{is}}{\sum_{i=1}^{n} e_{is}}
```

where $`e_{is}`$ is the expression level of gene $`i`$ at a given sample
$`s`$ (e.g. a biological replicate for a developmental stage), and
$`{ts}_i`$ is the evolutionary age of gene $`i`$.

Like TAI, for developmental time-course data, $`{TSI}_s`$ can be grouped
for each stage (if replicate data exists) and compared across
developmental stages to obtain the overall TSI profile.  
For pairwise comparisons (see
[📈](https://drostlab.github.io/myTAI/articles/tai-stats.md)),
$`{TSI}_s`$ can be grouped as one of the conditions to be compared.  
For single cell data, (see
[📊](https://drostlab.github.io/myTAI/articles/phylo-expression-object.md)),
$`{TSI}_s`$ is the TSI for a given cell or cell type (if pseudo-bulked).

Thus, TSI represents the expression-weighted mean gene expression
specificity.

You can re-use the same workflow as for TDI that we have shown above.
The only differences are that (1) you need to give each gene the
expression specificity score and not divergence strata, and (2) you will
be using real data not the mock dataset!

For TSI, accurate inference of expression specificity is key.

According to [Kryuchkova-Mostacci & Robinson-Rechavi,
2017](https://doi.org/10.1093/bib/bbw008):

> Tau appears consistently to be the most robust method in our analyses.

Thus, we recommend the *tau* (or ) score as introduced by [Yanai et
al. 2005](https://doi.org/10.1093/bioinformatics/bti042) as a measure of
expression specificity, though other measures exist (e.g. those based on
information entropy). You can also decile the *tau* values (similar to
deciled dNdS for divergence-strata) to enable analogous comparisons to
TAI.

Regardless of the measures chosen for expression specificity/breadth,
the estimates improve with data quality (i.e. number of RNA-seq samples)
and quantity (i.e. diversity of the sampling such as multiple tissues
and stages, and low or no technical noise).

## Transcriptome {insert your idea} index

Using two examples (TDI and TSI), we went through different (non-TAI)
evolutionary and expression indices that can be used for evolutionary
transcriptomics studies. As long as the computation follows the form
$`{TxI}_s = \frac{\sum_{i=1}^{n} {x}_i \cdot e_{is}}{\sum_{i=1}^{n} e_{is}}`$,
where $`x`$ is a value (factor or numeric) assigned to each gene, you
can even invent new ‘indices’.

Looking forward to seeing what can be done!

## Summary

You might wonder - why so many transcriptome indices? TAI, TDI, TSI, TPI
and other such indices provide **complementary** views on the
transcriptome evolution patterns. For more details on their biological
interpretation, please check these studies for
[`TAI`](https://www.nature.com/articles/nature09632),
[`TDI`](https://www.nature.com/articles/nature11394),
[`TSI`](https://www.nature.com/articles/s41586-024-08059-8),
[`TPI`](https://doi.org/10.1093/molbev/msw044), etc.

If you are using different indices, it is important to check that the
gene-level values assigned (e.g. age, divergence-stratum) are not
strongly correlated. High correlation means that the resulting indices
do **not provide complementary** information (i.e. not reflect different
biological aspects).
