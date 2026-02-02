# Getting started

  
  

  
  

  
  

## 1. Install **myTAI**

``` r
# from CRAN
install.packages("myTAI", dependencies = TRUE)

# or, the developer version containing the newest features
devtools::install_github("drostlab/myTAI")

# make sure myTAI version is > 2.0.0
packageVersion("myTAI")
```

  
  

  
  

## 2. Run **myTAI**

``` r
library(myTAI)
# obtain an example phylo-expression object
data("example_phyex_set")
# plot away!
myTAI::plot_signature(example_phyex_set)  
```

![plot_signature function
output](myTAI_files/figure-html/unnamed-chunk-2-1.png)

or with [your `BulkPhyloExpressionSet`
dataset](https://drostlab.github.io/myTAI/articles/phylo-expression-object.md).

![PhyEx overview](Figures/myTAI_phyex.png)

  
  

  
  

  
  

## 3. Enjoy 🍹

(and don’t forget to
[cite](https://doi.org/10.1093/bioinformatics/btx835))

  
  

  
  

  
  

  
  

Want to use myTAI on your data?  
→ click on the icons 🧚 below.

------------------------------------------------------------------------

[📊](https://drostlab.github.io/myTAI/articles/phylo-expression-object.md)
Bring your datasets into myTAI.  
[📈](https://drostlab.github.io/myTAI/articles/tai-stats.md) Statistical
analyses → learn about our permutation tests.  
[🛡️](https://drostlab.github.io/myTAI/articles/tai-transform.md) Check
the robustness of your patterns.  
[🔨](https://drostlab.github.io/myTAI/articles/tai-breaker.md) Destroy
the hourglass! → `gaTAI` functions to extract the genes that drive your
evolutionary transcriptomic patterns.  
[📚](https://drostlab.github.io/myTAI/articles/phylostratigraphy.md)
Read about gene age inference (phylostratigraphy) and transcriptome age
index.  
[🌄](https://drostlab.github.io/myTAI/articles/tai-gallery.md) Check out
our gallery of example plots + functions.

------------------------------------------------------------------------

  
  

  
  

  
  

Want to learn even more?  
→ check out our vignettes under the
[`Articles`](https://drostlab.github.io/myTAI/articles/index.md) tab or
the `Reference` tab for a list of functions and their documentation.

  
  

  
  

  
  
