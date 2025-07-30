#' Example phyex set
#'
#' Arabidopsis thaliana embryogenesis dataset from Hoffman et al. 2019. check the phyexSet package for more details
#'
#' @format A BulkPhyloExpressionSet with 8 developmental stages (3 reps each) and 27520 genes
#' @source phyexSets::Athaliana.embryogenesis_2019 matched with the Arabidopsis thaliana phylomap from phylomapr
"example_phyex_set"

#' Example phyex set old
#'
#' Arabidopsis thaliana embryogenesis dataset from Xiang et al. 2011. check the phyexSet package for more details
#'
#' @format A BulkPhyloExpressionObject with 7 developmental stages (1 rep each) and 25096 genes
#' @source phyexSets::Athaliana.embryogenesis_2011 matched with the Arabidopsis thaliana phylomap from phylomapr
"example_phyex_set_old"

#' Load Example Single-Cell PhyloExpressionSet
#'
#' Creates and returns an example ScPhyloExpressionSet object for testing and examples.
#' Only works if Seurat is installed.
#'
#' @return A ScPhyloExpressionSet object, or NULL if Seurat is not available.
#' @export
load_example_phyex_set_sc <- function() {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        warning("Seurat package is not installed. Returning NULL.")
        return(NULL)
    }

    set.seed(1234)
    n_genes <- 1000
    n_cells <- 1000

    counts <- matrix(
        stats::rnbinom(n_genes * n_cells, size = 5, mu = 100),
        nrow = n_genes,
        ncol = n_cells
    )
    rownames(counts) <- paste0("Gene-", 1:n_genes)
    colnames(counts) <- paste0("Cell-", 1:n_cells)

    metadata <- data.frame(
        groups = factor(sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE), 
                        levels = c("TypeA", "TypeB", "TypeC"), ordered = TRUE),
        nCount_RNA = colSums(counts),
        nFeature_RNA = colSums(counts > 0),
        row.names = colnames(counts)
    )

    example_seurat <- Seurat::CreateSeuratObject(
        counts = counts,
        meta.data = metadata,
        project = "ExampleSC"
    )
    Seurat::Idents(example_seurat) <- "groups"

    suppressWarnings({
        example_seurat <- Seurat::NormalizeData(example_seurat, verbose = FALSE)
        example_seurat <- Seurat::FindVariableFeatures(example_seurat, nfeatures = 50, verbose = FALSE)
        example_seurat <- Seurat::ScaleData(example_seurat, verbose = FALSE)
        example_seurat <- Seurat::RunPCA(example_seurat, npcs = 10, verbose = FALSE)
        example_seurat <- Seurat::RunUMAP(example_seurat, dims = 1:10, verbose = FALSE)
    })

    phylomap_example <- data.frame(
        Stratum = sample(1:10, nrow(example_seurat), replace = TRUE),
        GeneID = rownames(example_seurat)
    )

    example_phyex_set_sc <- match_map_sc(
        seurat = example_seurat,
        phylomap = phylomap_example,
        layer = "counts",
        name = "Single Cell Example"
    )

    return(example_phyex_set_sc)
}
