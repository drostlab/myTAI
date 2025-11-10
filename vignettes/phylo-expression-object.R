## ----message = FALSE, echo = FALSE--------------------------------------------
library(myTAI)
data("example_phyex_set")

## ----message = FALSE, eval=FALSE----------------------------------------------
# # inspect the object
# example_phyex_set

## ----message = FALSE, eval=FALSE----------------------------------------------
# S7::prop_names(example_phyex_set)

## ----message = FALSE, echo=FALSE----------------------------------------------
S7::prop_names(example_phyex_set)

## ----message = FALSE, eval = FALSE--------------------------------------------
# # let's explore the properties using @
# # for example:
# example_phyex_set@strata |> head()
# example_phyex_set@gene_ids |> head()
# example_phyex_set@expression[1:4,1:5]

## ----message = FALSE, eval = FALSE--------------------------------------------
# myTAI::stat_flatline_test(example_phyex_set)
# myTAI::stat_flatline_test(example_phyex_set)

## ----message = FALSE, echo=FALSE----------------------------------------------
data("example_phyex_set_old")
example_expression <- 
    example_phyex_set_old@expression |> 
    as.data.frame() |> 
    tibble::rownames_to_column(var = "GeneID")
example_phylorank <-
    example_phyex_set_old@strata

## ----message = FALSE, eval = FALSE--------------------------------------------
# # let's explore the example expression dataset
# # for example:
# example_expression |> head()

## ----message = FALSE, echo=FALSE----------------------------------------------
example_expression |> head()

## ----message = FALSE, echo = FALSE--------------------------------------------
example_phyex_set.df <-
    data.frame(phylorank = example_phylorank) |>
    tibble::rownames_to_column(var = "GeneID") |>
    dplyr::select(phylorank, GeneID) |>
    dplyr::left_join(example_expression, by = "GeneID")

## ----message = FALSE, eval=FALSE----------------------------------------------
# example_phyex_set.df |> head()
# example_phyex_set.df |> str()

## ----message = FALSE, echo=FALSE----------------------------------------------
example_phyex_set.df |> head()
example_phyex_set.df |> str()

## ----warning=FALSE, message=FALSE, echo = FALSE-------------------------------
example_phyex_set.remake <- 
    myTAI::BulkPhyloExpressionSet_from_df(
        example_phyex_set.df,
        name = "reconstituted phyex_set_old")

## ----message = FALSE, eval=FALSE----------------------------------------------
# example_phyex_set.remake

## ----message = FALSE, echo=FALSE----------------------------------------------
print(example_phyex_set.remake)

## ----message = FALSE, fig.height=4, fig.width=6, fig.alt="plot signature function output for bulk RNAseq", dev.args = list(bg = 'transparent'), fig.align='center'----
example_phyex_set.remake |> 
    myTAI::plot_signature()

## ----message = FALSE, echo = FALSE--------------------------------------------
data("example_phyex_set") # not data("example_phyex_set_old")
example_expression <- 
    example_phyex_set@expression |> 
    as.data.frame() |> 
    tibble::rownames_to_column(var = "GeneID")
example_phylorank <-
    example_phyex_set@strata

example_phyex_set.df <-
    data.frame(phylorank = example_phylorank) |>
    tibble::rownames_to_column(var = "GeneID") |>
    dplyr::select(phylorank, GeneID) |>
    dplyr::left_join(example_expression, by = "GeneID")

## ----message = FALSE, eval=FALSE----------------------------------------------
# colnames(example_phyex_set.df)

## ----message = FALSE, echo=FALSE----------------------------------------------
colnames(example_phyex_set.df)

## ----message = FALSE, echo = FALSE--------------------------------------------
groups_example_phyex_set.df <- c(
            rep("Preglobular",3), 
            rep("Globular",3), 
            rep("Early Heart",3), 
            rep("Late Heart",3), 
            rep("Early Torpedo",3), 
            rep("Late Torpedo",3), 
            rep("Bent Cotyledon",3), 
            rep("Mature Green",3))

# we can add the group information using groups
example_phyex_set.remake <- 
    myTAI::BulkPhyloExpressionSet_from_df(
        example_phyex_set.df,
        groups = groups_example_phyex_set.df, # adding group information here
        name = "reconstituted phyex_set")

## ----message = FALSE, warning=FALSE, eval = requireNamespace("Seurat", quietly = TRUE)----
# Load the relevant packages
library(dplyr)
library(Seurat)

pbmc_raw <- read.table(
  file = system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'),
  as.is = TRUE
)

## ----message = FALSE, eval=FALSE----------------------------------------------
# head(pbmc_raw)

## ----message = FALSE, eval = requireNamespace("Seurat", quietly = TRUE), echo=FALSE----
head(pbmc_raw)

## ----message = FALSE, warning = FALSE, eval = requireNamespace("Seurat", quietly = TRUE)----
pbmc_small <- CreateSeuratObject(counts = pbmc_raw)

## ----message = FALSE, eval = requireNamespace("Seurat", quietly = TRUE)-------
is(pbmc_small)

## ----message = FALSE, eval = requireNamespace("Seurat", quietly = TRUE)-------
gene_names <- pbmc_small |> rownames()
example_phylorank_sc <- 
    setNames(sample(1:10, length(gene_names), replace = TRUE), gene_names) |>
    as.factor()

## ----message = FALSE, eval = requireNamespace("Seurat", quietly = TRUE)-------
example_phyex_set_sc <- 
    myTAI::ScPhyloExpressionSet_from_seurat(
        pbmc_small,
        strata = example_phylorank_sc)

## ----message = FALSE, warning=FALSE, eval = requireNamespace("Seurat", quietly = TRUE)----
# example workflow to cluster the single cell data
pbmc_small <- Seurat::NormalizeData(pbmc_small)
pbmc_small <- Seurat::FindVariableFeatures(pbmc_small, selection.method = "vst", nfeatures = 20)
pbmc_small <- Seurat::ScaleData(pbmc_small)
pbmc_small <- Seurat::RunPCA(pbmc_small, features = VariableFeatures(object = pbmc_small))
pbmc_small <- Seurat::FindNeighbors(pbmc_small, dims = 1:10)
pbmc_small.cluster <- Seurat::FindClusters(pbmc_small, resolution = 0.8)

## ----message = FALSE, eval = requireNamespace("Seurat", quietly = TRUE)-------
example_phyex_set_sc.cluster <- 
    myTAI::ScPhyloExpressionSet_from_seurat(
        seurat = pbmc_small.cluster,
        strata = example_phylorank_sc)

## ----message = FALSE, eval = FALSE--------------------------------------------
# example_phyex_set_sc.cluster

## ----message = FALSE, fig.height=4, fig.width=6, fig.alt="plot signature function output for single cell", dev.args = list(bg = 'transparent'), fig.align='center', eval = requireNamespace("Seurat", quietly = TRUE)----
myTAI::plot_signature(example_phyex_set_sc.cluster)

