
#' @import dplyr
#' @export
get_sc_TAI <- function(seurat, phylomap, slot = "data") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }
    colnames(phylomap) <- c("Stratum", "GeneID")
    counts <- Seurat::GetAssayData(seurat, slot = slot) 
    gene_ids <- rownames(counts)
    phylo_map_sorted <- phylomap[match(gene_ids, phylomap$GeneID), ]
    ps_vec <- phylo_map_sorted |> pull(Stratum)
    ps_vec[is.na(ps_vec)] <- 0
    TAI <- (Matrix::t(counts) %*% ps_vec) / Matrix::colSums(counts)
    return(as.vector(TAI))
}