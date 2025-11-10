## ----message = FALSE, warning = FALSE-----------------------------------------
library(myTAI)

## ----message = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature function output no transform", dev.args = list(bg = 'transparent'), fig.align='center'----
data("example_phyex_set")
# no transformation
myTAI::plot_signature(example_phyex_set)

## ----message = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature function output rank transform", dev.args = list(bg = 'transparent'), fig.align='center'----
# rank transformation
myTAI::plot_signature(example_phyex_set |> myTAI::tf(FUN = function(x) apply(x, 2, base::rank)))

## ----message = FALSE, fig.height=5, fig.width=8, fig.alt="plot_signature_transformed function output transform all", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_signature_transformed(example_phyex_set)

## ----message = FALSE, warning = FALSE, results = FALSE------------------------
tf_stability_res <- myTAI::tf_stability(example_phyex_set)

## ----message = FALSE, warning = FALSE-----------------------------------------
tf_stability_res

## ----eval = FALSE, warning = FALSE--------------------------------------------
# COUNT_TRANSFORMS_REMOVED <- COUNT_TRANSFORMS[!names(COUNT_TRANSFORMS) %in% c("rlog", "rank", "vst")]

## ----eval = FALSE, warning = FALSE--------------------------------------------
# COUNT_TRANSFORMS_ADDED <- c(COUNT_TRANSFORMS, list(new_transform = function(x) x + 1))

## ----eval = FALSE, warning = FALSE--------------------------------------------
# myTAI::tf_stability(example_phyex_set, transformations = COUNT_TRANSFORMS_REMOVED)
# myTAI::plot_signature_transformed(example_phyex_set, transformations = COUNT_TRANSFORMS_REMOVED)

## ----message = FALSE, fig.height=5, fig.width=8, fig.alt="plot_signature_gene_quantiles function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_signature_gene_quantiles(example_phyex_set)

## ----message = FALSE----------------------------------------------------------
# get top 1% variance genes (default parameters)
genes.top_var <- myTAI::genes_top_variance(example_phyex_set)
# get top 1% expressed genes (default parameters)
genes.top_mean <- myTAI::genes_top_mean(example_phyex_set)
# get lowly expressed genes with mean < 1 (default parameters)
genes.low_expr <- myTAI::genes_lowly_expressed(example_phyex_set)

## ----message = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature function output removing top variance", dev.args = list(bg = 'transparent'), fig.align='center'----
# removing top 1% variance genes
example_phyex_set.rm_top_var <-
  myTAI::remove_genes(example_phyex_set, genes.top_var)
# plot
myTAI::plot_signature(example_phyex_set.rm_top_var)

## ----message = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature function output for top expressed", dev.args = list(bg = 'transparent'), fig.align='center'----
# select top 1% expressed genes
example_phyex_set.select_top_expr <- 
  myTAI::select_genes(example_phyex_set, genes.top_mean)
# plot
myTAI::plot_signature(example_phyex_set.select_top_expr)

## ----message = FALSE, fig.height=4, fig.width=6, fig.alt="plot_signature function output for gene of interest", dev.args = list(bg = 'transparent'), fig.align='center'----
# remove two genes (AT1G02780 & AT1G03880)
example_phyex_set.select_goi <- 
  myTAI::remove_genes(example_phyex_set, c("AT1G02780", "AT1G03880"))
# plot
myTAI::plot_signature(example_phyex_set.select_goi)

