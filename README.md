myTAI
=====

[![Travis-CI Build Status](https://travis-ci.org/HajkD/myTAI.svg?branch=master)](https://travis-ci.org/HajkD/myTAI)

### Performing Phylotranscriptomics with R

The `myTAI` package allows users to capture evolutionary signals in developmental transcriptomes using a phylotranscriptomic approach.
    
Phylotranscriptomics defines the concept of combining genetic sequence information with 
gene expression levels to capture evolutionary signals during development (Domazet-Loso and Tautz, 2010; Quint et al., 2012; Drost et al., 2015).

This subfield of Evolutionary Developmental Biology aims to investigate stages or periods
of evolutionary conservation or constraint in developmental processes of extant species.

`myTAI` provides easy to use and optimized functions to perform phylostrancriptomic analyses for any annotated organism and developmental process of interest. Additionally, a customized visualization framework implemented in myTAI allows users to generate publication quality plots of their custom phylotranscriptomic analyses.


## Tutorials

These tutorials will get you started with this package:

- [A Brief Introduction to the myTAI Package](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd)
- [Intermediate Concepts of Phylotranscriptomics](https://github.com/HajkD/myTAI/blob/master/vignettes/Intermediate.Rmd)
- [Advanced Topics of Phylotranscriptomics](https://github.com/HajkD/myTAI/blob/master/vignettes/Advanced.Rmd)
- [Perform Age Enrichment Analyses](https://github.com/HajkD/myTAI/blob/master/vignettes/Enrichment.Rmd)
- [Gene Expression Analysis with myTAI](https://github.com/HajkD/myTAI/blob/master/vignettes/Expression.Rmd)
- [Taxonomic Information Retrieval](https://github.com/HajkD/myTAI/blob/master/vignettes/Taxonomy.Rmd)

You can also read the tutorials within ([RStudio](http://www.rstudio.com/)) :

```r

# first install the myTAI package 
# -> see "Fast Installation Guide" for the current development version
install.packages("myTAI", repos = "https://cran.rstudio.com/", dependencies = TRUE, type = "source")

# to perform differential gene expression analyses with myTAI
# please install the edgeR package
# install edgeR
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

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

In the `myTAI` framework you can find:

#### Phylotranscriptiomics Measures:

* `TAI()` : Function to compute the Transcriptome Age Index (TAI)
* `TDI()` : Function to compute the Transcriptome Divergence Index (TDI)
* `REMatrix()` : Function to compute the relative expression profiles of all phylostrata or divergence-strata
* `RE()` : Function to transform mean expression levels to relative expression levels
* `pTAI()` : Compute the Phylostratum Contribution to the global TAI
* `pTDI()` : Compute the Divergence Stratum Contribution to the global TDI
* `pMatrix()` : Compute Partial TAI or TDI Values
* `pStrata()` : Compute Partial Strata Values

#### Visualization and Analytics Tools:

* `PlotPattern()` : Function to plot the TAI or TDI profiles and perform statistical tests
* `PlotCorrelation()` : Function to plot the correlation between phylostratum values and divergence-stratum values
* `PlotRE()` : Function to plot the relative expression profiles
* `PlotBarRE()` : Function to plot the mean relative expression levels of phylostratum or divergence-stratum classes as barplot
* `PlotMeans()` : Function to plot the mean expression profiles of phylostrata or divergence-strata
* `PlotDistribution()` : Function to plot the frequency distribution of genes within the corresponding phylostratigraphic map or divergence map
* `PlotContribution()` : Plot the Phylostratum or Divergence Stratum Contribution to the Global TAI/TDI Pattern
* `PlotEnrichment()` : Plot the Phylostratum or Divergence Stratum Enrichment of a given Gene Set
* `PlotGeneSet()` : Plot the Expression Profiles of a Gene Set


#### A Statistical Framework and Test Statistics:

* `FlatLineTest()` : Function to perform the __Flat Line Test__ that quantifies the statistical significance of an observed
phylotranscriptomics pattern (significant deviation from a frat line = evolutionary signal)
* `ReductiveHourglassTest()` : Function to perform the __Reductive Hourglass Test__ that statistically evaluates the existence of a phylotranscriptomic hourglass pattern (hourglass model)
* `EarlyConservationTest()` : Function to perform the __Reductive Early Conservation Test__ that statistically evaluates the existence of a monotonically increasing phylotranscriptomic pattern (early conservation model)
* `EnrichmentTest()` : Phylostratum or Divergence Stratum Enrichment of a given Gene Set based on Fisher's Test
* `bootMatrix()` : Compute a Permutation Matrix for Test Statistics

All three functions also include visual analytics tools to quantify the goodness of test statistics.

#### Differential Gene Expression Analysis

* `DiffGenes()` : Implements Popular Methods for Differential Gene Expression Analysis
* `CollapseReplicates()` : Combine Replicates in an ExpressionSet
* `CombinatorialSignificance()` : Compute the Statistical Significance of Each Replicate Combination
* `Expressed()` : Filter Expression Levels in Gene Expression Matrices (define expressed genes)
* `SelectGeneSet()` : Select a Subset of Genes in an ExpressionSet

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


## Install Developer Version of myTAI

```r
# The developer version can be installed directly from github:

# install.packages("devtools")

# install developer version of myTAI
library(devtools)
install_github("HajkD/myTAI", build_vignettes = TRUE, dependencies = TRUE)

# On Windows, this won't work - see ?build_github_devtools
# install_github("HajkD/myTAI", build_vignettes = TRUE, dependencies = TRUE)

# When working with Windows, first you need to install the
# R package: rtools -> http://cran.r-project.org/bin/windows/Rtools/
# or consult: http://www.rstudio.com/products/rpackages/devtools/

# Afterwards you can install devtools -> install.packages("devtools")
# and then you can run:

devtools::install_github("HajkD/myTAI", build_vignettes = TRUE, dependencies = TRUE)

# and then call it from the library
library("myTAI", lib.loc = "C:/Program Files/R/R-3.1.1/library")

```

## References

Domazet-Lošo T. and Tautz D. __A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns__. _Nature_ (2010) 468: 815-8.

Quint M. et al. __A transcriptomic hourglass in plant embryogenesis__. _Nature_ (2012) 490: 98-101.

Drost HG, Gabel A, Grosse I, Quint M. __Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis__. _Mol. Biol. Evol._ (2015) 32 (5): 1221-1231.

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


