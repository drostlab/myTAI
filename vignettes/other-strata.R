## ----message = FALSE----------------------------------------------------------
library(myTAI)
data("example_phyex_set_old")

## ----message = FALSE----------------------------------------------------------
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

## ----message = FALSE----------------------------------------------------------
example_phyex_set_tdi.df <-
    data.frame(phylorank = example_divergence_strata) |>
    tibble::rownames_to_column(var = "GeneID") |>
    dplyr::select(phylorank, GeneID) |>
    dplyr::left_join(example_expression, by = "GeneID")

## ----message = FALSE----------------------------------------------------------
example_phyex_set_tdi <- 
    myTAI::as_BulkPhyloExpressionSet(
        example_phyex_set_tdi.df,
        name = "Example TDI data",
        index_type = "TDI")

## ----message = FALSE, fig.height=4, fig.width=6, fig.alt="plot signature function output for mock TDI data", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::plot_signature(example_phyex_set_tdi)

## ----message = FALSE, warning = FALSE, results = FALSE, fig.height=3, fig.width=5, fig.alt="stat_flatline_test function output", dev.args = list(bg = 'transparent'), fig.align='center'----
myTAI::stat_flatline_test(example_phyex_set_tdi)

