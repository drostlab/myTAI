# Create example single-cell PhyloExpressionSet for testing

# Load required libraries

set.seed(1234)

if (requireNamespace("Seurat", quietly = TRUE)) {
    library(Seurat)
    # Create a synthetic Seurat object for testing
    n_genes <- 100
    n_cells <- 200
    
    # Create synthetic expression data with some realistic patterns
    counts <- matrix(
        rpois(n_genes * n_cells, lambda = 5),
        nrow = n_genes,
        ncol = n_cells
    )
    rownames(counts) <- paste0("Gene-", 1:n_genes)
    colnames(counts) <- paste0("Cell-", 1:n_cells)
    
    # Create synthetic metadata
    metadata <- data.frame(
        groups = factor(sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE)),
        nCount_RNA = colSums(counts),
        nFeature_RNA = colSums(counts > 0),
        row.names = colnames(counts)
    )
    
    # Create Seurat object
    example_seurat <- Seurat::CreateSeuratObject(
        counts = counts,
        meta.data = metadata,
        project = "ExampleSC"
    )
    
    # Process the data: normalize, find variable features, scale, PCA, and UMAP
    suppressWarnings({
        example_seurat <- Seurat::NormalizeData(example_seurat, verbose = FALSE)
        example_seurat <- Seurat::FindVariableFeatures(example_seurat, nfeatures = 50, verbose = FALSE)
        example_seurat <- Seurat::ScaleData(example_seurat, verbose = FALSE)
        example_seurat <- Seurat::RunPCA(example_seurat, npcs = 10, verbose = FALSE)
        example_seurat <- Seurat::RunUMAP(example_seurat, dims = 1:10, verbose = FALSE)
    })
    
    # Create a dummy phylomap for testing
    phylomap_example <- data.frame(
        Stratum = sample(1:10, nrow(example_seurat), replace = TRUE),
        GeneID = rownames(example_seurat)
    )
    
    example_phyex_set_sc <- as_ScPhyloExpressionSet(
        seurat = example_seurat,
        phylomap = phylomap_example,
        cell_identity = "groups",
        name = "Single Cell Example"
    )
}

# Save the data
usethis::use_data(example_phyex_set_sc, overwrite = TRUE)
