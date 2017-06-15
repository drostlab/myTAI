myTAI
=====
[![Travis-CI Build Status](https://travis-ci.org/HajkD/myTAI.svg?branch=master)](https://travis-ci.org/HajkD/myTAI)
 [![rpackages.io rank](https://www.rpackages.io/badge/myTAI.svg)](https://www.rpackages.io/package/myTAI)

### Evolutionary Transcriptomics with R

Today, phenotypic phenomena such as morphological mutations, diseases or developmental processes are primarily investigated on the molecular level using transcriptomics approaches. Transcriptomes denote the total number of quantifiable transcripts present at a specific stage in a biological process. In disease or developmental (defect) studies transcriptomes are usually measured over several time points. In treatment studies aiming to quantify differences in the transcriptome due to biotic stimuli, abiotic stimuli, or diseases usually treatment / disease versus non-treatment / non-disease transcriptomes are being compared. In either case, comparing changes in transcriptomes over time or between treatments allows us to identify genes and gene regulatory mechanisms that might be involved in governing the biological process of investigation. Although transcriptomics studies are based on a powerful methodology little is known about the evolution of such transcriptomes. Understanding the evolutionary mechanism that change transcriptomes over time, however, might give us a new perspective on how diseases emerge in the first place or how morphological changes are triggered by changes of developmental transcriptomes.

Evolutionary transcriptomics aims to capture and quantify the evolutionary conservation of genes that contribute to the transcriptome during a specific stage of the biological process of interest. This quantification on the highest level is achieved through transcriptome indices ([Domazet-Lošo and Tautz, 2010](http://www.nature.com/nature/journal/v468/n7325/abs/nature09632.html); [Drost et al., 2016a](http://biorxiv.org/content/early/2016/05/03/051565)) which denote weighted means of gene age or rate of protein substitutions. In general, evolutionary transcriptomics can be used as a method to quantify the evolutionary conservation of transcriptomes to investigate how transcriptomes underlying biological processes are constrained or channeled due to evolutionary history (Dollow's law) ([Drost et al., 2017](http://www.sciencedirect.com/science/article/pii/S0959437X16302040)).

In principle, any transcriptome dataset published so far can be combined with evolutionary information. Thus, `myTAI` in combination with evolutionary information can be used to study corresponding transcriptomes with any available transcriptome dataset. 

For the purpose of performing large scale evolutionary transcriptomics studies, the `myTAI` package implements frameworks to allow researchers to study the evolution of biological processes and to detect stages or periods of evolutionary conservation or variability. 

I hope that `myTAI` will become the community standard tool to perform evolutionary transcriptomics studies and I am happy to add required functionality upon request.


The following tutorials will provide use cases and detailed explainations of how to quantify transcriptome onservation with `myTAI` and how to interpret the results generated with this software tool.

### Citation

**Please cite one of the following references when using `myTAI` for your own research. This will allow me to continue
working on this software tool and will motivate me to extend its functionality and usability. Many thanks in advance :)**

>  Drost HG, Gabel A, Domazet-Lošo T, Grosse I, Quint M. 2016. __Capturing Evolutionary Signatures in Transcriptomes with myTAI__. 
[doi: https://doi.org/10.1101/051565](http://biorxiv.org/content/early/2016/05/03/051565)
>
> Drost HG, Gabel A, Grosse I, Quint M. 2015. __Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis__. _Mol. Biol. Evol._ 32 (5): 1221-1231. [doi:10.1093/molbev/msv012](http://mbe.oxfordjournals.org/content/32/5/1221.abstract?sid=767aea12-1eb3-40be-8c6a-e2861f159b46)


## Installation

Users can download `myTAI` from CRAN :

```r
# install myTAI 0.5.0 from CRAN
source("http://bioconductor.org/biocLite.R")
biocLite('myTAI')
```

## Install Developer Version
Some bug fixes or new functionality will not be available on CRAN yet, but in
the developer version here on GitHub. To download and install the most recent
version of `myTAI` run:

```r
# install the developer version of myTAI on your system
source("http://bioconductor.org/biocLite.R")
biocLite("HajkD/myTAI")
```

## NEWS

The current status of the package as well as a detailed history of the
functionality of each version of `myTAI` can be found in the [NEWS](https://github.com/HajkD/myTAI/blob/master/NEWS.md) section.



## Tutorials

These tutorials introduce users to `myTAI`:

- [Introduction to the myTAI Package](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd)
- [Intermediate Concepts of Phylotranscriptomics](https://github.com/HajkD/myTAI/blob/master/vignettes/Intermediate.Rmd)
- [Advanced Topics of Phylotranscriptomics](https://github.com/HajkD/myTAI/blob/master/vignettes/Advanced.Rmd)
- [Perform Age Enrichment Analyses](https://github.com/HajkD/myTAI/blob/master/vignettes/Enrichment.Rmd)
- [Gene Expression Analysis with myTAI](https://github.com/HajkD/myTAI/blob/master/vignettes/Expression.Rmd)
- [Taxonomic Information Retrieval](https://github.com/HajkD/myTAI/blob/master/vignettes/Taxonomy.Rmd)


### Package Dependencies


```r
# to perform differential gene expression analyses with myTAI
# please install the edgeR package
# install edgeR
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
```

## Getting started with `myTAI`

Users can also read the tutorials within ([RStudio](http://www.rstudio.com/)) :

```r
# source the myTAI package
library(myTAI)

# look for all tutorials (vignettes) available in the myTAI package
# this will open your web browser
browseVignettes("myTAI")

# or as single tutorials

# open tutorial: Introduction to Phylotranscriptomics and myTAI
 vignette("Introduction", package = "myTAI")

# open tutorial: Intermediate Concepts of Phylotranscriptomics
 vignette("Intermediate", package = "myTAI")

# open tutorial: Advanced Concepts of Phylotranscriptomics
 vignette("Advanced", package = "myTAI")

# open tutorial: Age Enrichment Analyses
 vignette("Enrichment", package = "myTAI")
 
# open tutorial: Gene Expression Analysis with myTAI
 vignette("Expression", package = "myTAI")
 
 # open tutorial: Taxonomic Information Retrieval with myTAI
 vignette("Taxonomy", package = "myTAI")
```

In the `myTAI` framework users can find:

#### Phylotranscriptomics Measures:

* `TAI()` : Function to compute the Transcriptome Age Index (TAI)
* `TDI()` : Function to compute the Transcriptome Divergence Index (TDI)
* `TPI()` : Function to compute the Transcriptome Polymorphism Index (TPI)
* `REMatrix()` : Function to compute the relative expression profiles of all phylostrata or divergence-strata
* `RE()` : Function to transform mean expression levels to relative expression levels
* `pTAI()` : Compute the Phylostratum Contribution to the global TAI
* `pTDI()` : Compute the Divergence Stratum Contribution to the global TDI
* `pMatrix()` : Compute Partial TAI or TDI Values
* `pStrata()` : Compute Partial Strata Values

#### Visualization and Analytics Tools:

* `PlotSignature()` : Main visualization function to plot evolutionary signatures across transcriptomes
* `PlotPattern()` : Base graphics function to plot evolutionary signatures across transcriptomes
* `PlotContribution()` : Plot Cumuative Transcriptome Index
* `PlotCorrelation()` : Function to plot the correlation between phylostratum values and divergence-stratum values
* `PlotRE()` : Function to plot the relative expression profiles
* `PlotBarRE()` : Function to plot the mean relative expression levels of phylostratum or divergence-stratum classes as barplot
* `PlotMeans()` : Function to plot the mean expression profiles of age categories
* `PlotMedians()` : Function to plot the median expression profiles of age categories
* `PlotVars()` : Function to plot the expression variance profiles of age categories
* `PlotDistribution()` : Function to plot the frequency distribution of genes within the corresponding age categories
* `PlotCategoryExpr()` : Plot the Expression Levels of each Age or Divergence Category as Barplot or Violinplot
* `PlotEnrichment()` : Plot the Phylostratum or Divergence Stratum Enrichment of a given Gene Set
* `PlotGeneSet()` : Plot the Expression Profiles of a Gene Set
* `PlotGroupDiffs()` : Plot the significant differences between gene expression distributions of PS or DS groups
* `PlotSelectedAgeDistr()` : Plot the PS or DS distribution of a selected set of genes

#### A Statistical Framework and Test Statistics:

* `FlatLineTest()` : Function to perform the __Flat Line Test__ that quantifies the statistical significance of an observed
phylotranscriptomics pattern (significant deviation from a frat line = evolutionary signal)
* `ReductiveHourglassTest()` : Function to perform the __Reductive Hourglass Test__ that statistically evaluates the existence of a phylotranscriptomic hourglass pattern (hourglass model)
* `EarlyConservationTest()` : Function to perform the __Reductive Early Conservation Test__ that statistically evaluates the existence of a monotonically increasing phylotranscriptomic pattern (early conservation model)
* `EnrichmentTest()` : Phylostratum or Divergence Stratum Enrichment of a given Gene Set based on Fisher's Test
* `bootMatrix()` : Compute a Permutation Matrix for Test Statistics

All functions also include visual analytics tools to quantify the goodness of test statistics.

#### Differential Gene Expression Analysis

* `DiffGenes()` : Implements Popular Methods for Differential Gene Expression Analysis
* `CollapseReplicates()` : Combine Replicates in an ExpressionSet
* `CombinatorialSignificance()` : Compute the Statistical Significance of Each Replicate Combination
* `Expressed()` : Filter Expression Levels in Gene Expression Matrices (define expressed genes)
* `SelectGeneSet()` : Select a Subset of Genes in an ExpressionSet
* `PlotReplicateQuality()` : Plot the Quality of Biological Replicates
* `GroupDiffs()` : Quantify the significant differences between gene expression distributions of PS or DS groups

#### Taxonomic Information Retrieval

* `taxonomy()` : Retrieve Taxonomic Information for any Organism of Interest

#### Minor Functions for Better Usibility and Additional Analyses

* `MatchMap()` : Match a Phylostratigraphic Map or Divergence Map with a ExpressionMatrix
* `tf()` : Transform Gene Expression Levels
* `age.apply()` : Age Category Specific apply Function
* `ecScore()` : Compute the Hourglass Score for the EarlyConservationTest
* `geom.mean()` : Geometric Mean
* `harm.mean()` : Harmonic Mean
* `omitMatrix()` : Compute TAI or TDI Profiles Omitting a Given Gene
* `rhScore()` : Compute the Hourglass Score for the Reductive Hourglass Test


## Developer Version of `myTAI`

The developer version of `myTAI` might include more functionality than the stable version on CRAN.
Hence users can download the current developer version of `myTAI` by typing:

```r
# The developer version can be installed directly from github:

# install.packages("devtools")

# install developer version of myTAI
library(devtools)
install_github("HajkD/myTAI", build_vignettes = TRUE, dependencies = TRUE)

# On Windows, this won't work - see ?build_github_devtools
# install_github("HajkD/myTAI", build_vignettes = TRUE, dependencies = TRUE)

# When working with Windows, first you need to install the
# R package: rtools 
# or consult: http://www.rstudio.com/products/rpackages/devtools/

# Afterwards you can install devtools -> install.packages("devtools")
# and then you can run:

devtools::install_github("HajkD/myTAI", build_vignettes = TRUE, dependencies = TRUE)

# and then call it from the library
library("myTAI", lib.loc = "C:/Program Files/R/R-3.1.1/library")

```

## References

Domazet-Lošo T. and Tautz D. __A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns__. _Nature_ (2010) 468: 815-8.

Quint M, Drost HG, et al. __A transcriptomic hourglass in plant embryogenesis__. _Nature_ (2012) 490: 98-101.

Drost HG, Gabel A, Grosse I, Quint M. __Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis__. _Mol. Biol. Evol._ (2015) 32 (5): 1221-1231.

Drost HG, Bellstädt J, Ó'Maoiléidigh DS, Silva AT, Gabel A, Weinholdt C, Ryan PT, Dekkers BJW, Bentsink L, Hilhorst H, Ligterink W, Wellmer F, Grosse I, and Quint M. __Post-embryonic hourglass patterns mark ontogenetic transitions in plant development__. _Mol. Biol. Evol._ (2016) [doi:10.1093/molbev/msw039](http://mbe.oxfordjournals.org/content/early/2016/02/23/molbev.msw039.short?rss=1) 

## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

https://github.com/HajkD/myTAI/issues


## Acknowledgement

I would like to thank several individuals for making this project possible.

First I would like to thank Ivo Grosse and Marcel Quint for providing me a place
and the environment to be able to work on fascinating topics of Evo-Devo research and for the
fruitful discussions that led to projects like this one.

Furthermore, I would like to thank Alexander Gabel and Jan Grau for valuable discussions
on how to improve some methodological concepts of some analyses present in this package.

I would also like to thank Master Students: Sarah Scharfenberg, Anne Hoffmann, and Sebastian Wussow
who worked intensively with this package and helped me to improve the usability and logic of the package environment.


