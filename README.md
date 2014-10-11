myTAI
=====

### A package to perform phylotranscriptomics analyses and visualization for Evolutionary Developmental Biology research.

The present collection of [R](http://cran.r-project.org/) functions can be used to perform phylotranscriptomics 
analyses and visualization to investigate phenomena within the field of Evolutionary Developmental Biology.
    
Phylotranscriptomics defines the concept of combining genetic sequence information with 
gene expression levels to capture evolutionary signals during developmental processes.

## Fast installation guide

```r
# install.packages("devtools")

# install the current version of myTAI on your system
library(devtools)
install_github("HajkD/myTAI")

# On Windows, this won't work - see ?build_github_devtools
install_github("HajkD/myTAI")

# When working with Windows, first you need to install the
# R package: rtools -> install.packages("rtools")

# Afterwards you can install devtools -> install.packages("devtools")
# and then you can run:

devtools::install_github("HajkD/myTAI")

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


You can also read the tutorials within R ([RStudio](http://www.rstudio.com/)) :

```r

# first install the myTAI package -> see "Fast Installation Guide" for the current development version
install.packages("myTAI", dependencies = TRUE)

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



## Acknowledgement

I would like to thank several individuals for making this project possible.

First I would like to thank Ivo Grosse and Marcel Quint for providing me a place
and the environment to be able to work on fascinating topics of Evo-Devo research and for the
fruitful discussions that lead to projects like this one.

Furthermore, I would like to thank Alexander Gabel and Jan Grau for valuable discussions
on how to improve some methodological concepts of some analyses present in this package.

I would also like to thank the Master of Science Students: Sarah Scharfenberg, Anne Hoffmann, and Sebastian Wussow
who worked intensively with this package and helped me to improve the usability and logic of this package.

As google scholar properly cites: [_On the shoulders of giants_](http://scholar.google.com/), I couldn't agree more than to add that this work
wouldn't be possible without the _giants_ that came before me and who provided such amazing work I can now benefit from and
I hope others can benefit from too.



