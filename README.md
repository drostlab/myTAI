myTAI
=====

### A Framework to Perform Phylotranscriptomics Analyses and Visualization for Evolutionary Developmental Biology Research.

The present collection of [R](http://cran.r-project.org/) functions can be used to perform phylotranscriptomics 
analyses and visualization to investigate phenomena within the field of Evolutionary Developmental Biology.
    
Phylotranscriptomics defines the concept of combining genetic sequence information with 
gene expression levels to capture evolutionary signals during developmental processes ([Domazet-Loso and Tautz, 2010](http://www.nature.com/nature/journal/v468/n7325/full/nature09632.html); [Quint et al., 2012](http://www.nature.com/nature/journal/v490/n7418/full/nature11394.html)).


In the `myTAI` framework you can find:

#### The following phylotranscriptiomics measures:

* `TAI()` : Function to compute the Transcriptome Age Index (TAI)
* `TDI()` : Function to compute the Transcriptome Divergence Index (TDI)
* `REMatrix()` : Function to compute the relative expression profiles of all phylostrata or divergence-strata

#### The following visualization and analytics tools:

* `PlotPattern()` : Function to plot the TAI or TDI profiles and perform statistical tests
* `PlotCorrelation()` : Function to plot the correlation between phylostratum values and divergence-stratum values
* `PlotRE()` : Function to plot the relative expression profiles
* `PlotBarRE()` : Function to plot the mean relative expression levels of phylostratum or divergence-stratum classes as barplot
* `PlotMeans()` : Function to plot the mean expression profiles of phylostrata or divergence-strata
* `PlotDistribution()` : Function to plot the frequency distribution of genes within the corresponding phylostratigraphic map or divergence map


#### The following statistical framework and test statistics:

* `FlatLineTest()` : Function to perform the __Flat Line Test__ that quantifies the statistical significance of an observed
phylotranscriptomics pattern (significant deviation from a frat line = no evolutionary signal)
* `ReductiveHourglassTest()` : Function to perform the __Reductive Hourglass Test__ that statistically evaluates the existence of a phylotranscriptomic hourglass pattern (hourglass model)
* `EarlyConservationTest()` : Function to perform the __Reductive Early Conservation Test__ that statistically evaluates the existence of a monotonically increasing phylotranscriptomic pattern (early conservation model)

All three functions also include visual analytics tools to quantify the goodness of test statistics.

#### Minor functions for better usibility and additional analytics

`age.apply()`, `bootMatrix()`, `combinatorialSignificance()`, `ecScore()`, `MatchMap()`, `omitMatrix()`, `pMatrix()`,
`RE()`, `rhScore()`, and `tf()`


## Fast installation guide

```r

# install myTAI version 0.0.1 from CRAN
install.packages("myTAI")


# The developer version can be installed directly from github:

# install.packages("devtools")

# install the current version of myTAI on your system
# this can take some time since the vignettes are
# very comprehensive and take some time to build
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

This package provides a broad toolbox for common phylotranscriptomics analyses and also includes statistical tests to verify observed phenomena and transcriptomic patterns.

Using the [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html) package, 
all computationally expansive functions have been written in C++ 
to enable fast analytics on phylotranscriptomics datasets.


## Tutorials

Three tutorials will get you started with this package:

- [An Introduction to the myTAI package](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd)
- [Intermediate concepts of the myTAI package](https://github.com/HajkD/myTAI/blob/master/vignettes/Intermediate.Rmd)
- [Advanced topics using the myTAI package](https://github.com/HajkD/myTAI/blob/master/vignettes/Advanced.Rmd)


You can also read the tutorials within ([RStudio](http://www.rstudio.com/)) :

```r

# first install the myTAI package 
# -> see "Fast Installation Guide" for the current development version
install.packages("myTAI", build_vignettes = TRUE, dependencies = TRUE)

# source the myTAI package
library(myTAI)

# look for all tutorials (vignettes) available in the myTAI package
# this will open your web browser
browseVignettes("myTAI")

# or as single tutorials

# open tutorial: Introduction
 vignette("Introduction", package = "myTAI")

# open tutorial: Intermediate
 vignette("Intermediate", package = "myTAI")

# open tutorial: Advanced
 vignette("Advanced", package = "myTAI")


```


## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

hajk-georg.drost@informatik.uni-halle.de


## Acknowledgement

I would like to thank several individuals for making this project possible.

First I would like to thank Ivo Grosse and Marcel Quint for providing me a place
and the environment to be able to work on fascinating topics of Evo-Devo research and for the
fruitful discussions that led to projects like this one.

Furthermore, I would like to thank Alexander Gabel and Jan Grau for valuable discussions
on how to improve some methodological concepts of some analyses present in this package.

I would also like to thank the Master Students: Sarah Scharfenberg, Anne Hoffmann, and Sebastian Wussow
who worked intensively with this package and helped me to improve the usability and logic of the package environment.


