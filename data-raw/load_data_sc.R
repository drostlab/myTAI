library(myTAI)

load_example_phyex_set_sc <- function() {
    set.seed(1234)
    n_genes <- 1000
    n_cells <- 1000

    # Create count matrix
    counts <- matrix(
        stats::rnbinom(n_genes * n_cells, size = 5, mu = 100),
        nrow = n_genes,
        ncol = n_cells
    )
    rownames(counts) <- paste0("Gene-", 1:n_genes)
    colnames(counts) <- paste0("Cell-", 1:n_cells)

    # Create metadata
    metadata <- data.frame(
        groups = factor(sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE), 
                        levels = c("TypeA", "TypeB", "TypeC"), ordered = TRUE),
        day = factor(sample(c("Day1", "Day3", "Day5", "Day7"), n_cells, replace = TRUE),
                     levels = c("Day1", "Day3", "Day5", "Day7"), ordered = TRUE),
        condition = factor(sample(c("Control", "Treatment"), n_cells, replace = TRUE),
                          levels = c("Control", "Treatment")),
        batch = factor(sample(c("Batch1", "Batch2", "Batch3"), n_cells, replace = TRUE)),
        nCount_RNA = colSums(counts),
        nFeature_RNA = colSums(counts > 0),
        row.names = colnames(counts)
    )

    # Create phylomap
    phylomap_example <- data.frame(
        Stratum = sample(1:10, n_genes, replace = TRUE),
        GeneID = rownames(counts)
    )

    # Create ScPhyloExpressionSet using the matrix constructor
    strata_values <- phylomap_example$Stratum[match(rownames(counts), phylomap_example$GeneID)]
    strata_factor <- factor(strata_values, levels = sort(unique(strata_values)))
    names(strata_factor) <- rownames(counts)
    
    example_phyex_set_sc <- ScPhyloExpressionSet_from_matrix(
        expression_matrix = Matrix::Matrix(counts, sparse = TRUE),
        strata = strata_factor,
        metadata = metadata,
        groups_column = "groups",
        name = "Single Cell Example",
        species = "Example Species",
        index_type = "TXI",
        identities_label = "Cell Type"
    )

    return(example_phyex_set_sc)
}

example_phyex_set_sc <- load_example_phyex_set_sc()

usethis::use_data(example_phyex_set_sc, overwrite = TRUE)
