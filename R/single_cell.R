
#' @title Calculate Single-Cell Transcriptomic Age Index
#' @description Compute the Transcriptomic Age Index (TAI) for individual cells
#' in a Seurat object using phylostratum information.
#' 
#' @param seurat A Seurat object containing single-cell expression data
#' @param phylomap A data frame with phylostratum information (columns: Stratum, GeneID)
#' @param slot Character string specifying which slot to use from the Seurat object (default: "data")
#' 
#' @return A numeric vector of TAI values for each cell
#' 
#' @details
#' This function calculates the Transcriptomic Age Index for each cell by:
#' 1. Extracting expression data from the specified Seurat slot
#' 2. Mapping genes to their phylostratum assignments
#' 3. Computing weighted averages of phylostratum values for each cell
#' 
#' Genes without phylostratum information are assigned a value of 0.
#' 
#' @examples
#' # Calculate TAI for single cells
#' # seurat_obj$TAI <- get_sc_TAI(seurat_obj, phylo_map)
#' 
#' @import dplyr
#' @export
get_sc_TAI <- function(seurat, phylomap, slot = "data") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }
    colnames(phylomap) <- c("Stratum", "GeneID")
    
    # Handle both old and new Seurat API
    counts <- tryCatch({
        # Try new API first (layer parameter)
        if (slot == "data") {
            # Try to get data layer, fall back to counts if empty
            tryCatch({
                result <- Seurat::GetAssayData(seurat, layer = "data")
                if (nrow(result) == 0 || ncol(result) == 0) {
                    # Fall back to counts if data layer is empty
                    Seurat::GetAssayData(seurat, layer = "counts")
                } else {
                    result
                }
            }, error = function(e) {
                # If data layer doesn't exist, use counts
                Seurat::GetAssayData(seurat, layer = "counts")
            })
        } else {
            Seurat::GetAssayData(seurat, layer = slot)
        }
    }, error = function(e) {
        # Fall back to old API
        Seurat::GetAssayData(seurat, slot = slot)
    })
    
    gene_ids <- rownames(counts)
    phylo_map_sorted <- phylomap[match(gene_ids, phylomap$GeneID), ]
    ps_vec <- phylo_map_sorted |> pull(Stratum)
    ps_vec[is.na(ps_vec)] <- 0
    TAI <- (Matrix::t(counts) %*% ps_vec) / Matrix::colSums(counts)
    return(as.vector(TAI))
}