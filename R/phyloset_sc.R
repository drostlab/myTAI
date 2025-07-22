#' @title Single-Cell PhyloExpressionSet S7 Class
#' @description S7 class for storing and manipulating single-cell phylotranscriptomic expression data.
#' This class extends PhyloExpressionSet where each cell is treated as a replicate and cell types
#' are the conditions after pseudobulking.
#' 
#' @param seurat A Seurat object containing single-cell expression data
#' @param cell_identity Character string specifying which cell identity to use from Seurat metadata
#' @param slot Character string specifying which slot to use from the Seurat object (default: "data")
#' @param min_cells_per_type Minimum number of cells required per cell type (default: 10)
#' 
#' @details
#' The ScPhyloExpressionSet class inherits all properties from PhyloExpressionSet:
#' - `counts`: Raw single-cell expression matrix from Seurat object
#' - `counts_collapsed`: Pseudobulked expression data across cell types  
#' - `TXI_reps`: Cell-level TAI values (each cell as a replicate)
#' - `TXI`: Cell type-level TAI values (pseudobulked)
#' - `groups`: Factor mapping each cell to its cell type
#' - `conditions`: Factor of unique cell type names
#' 
#' Additional single-cell specific properties:
#' - `cell_metadata`: Cell-level metadata from Seurat object
#' 
#' @examples
#' # Create a ScPhyloExpressionSet from Seurat object
#' # sc_phyex_set <- as_ScPhyloExpressionSet(seurat_obj, phylo_map, "celltype")
#' 
#' @import S7
#' @export
ScPhyloExpressionSet <- new_class("ScPhyloExpressionSet",
    parent = PhyloExpressionSet,
    properties = list(
        seurat = new_required_property(
            name = "seurat"
        ),
        cell_identity = new_required_property(
            class = class_character,
            name = "cell_identity"
        ),
        slot = new_property(
            class = class_character,
            default = "data"
        ),
        min_cells_per_type = new_property(
            class = class_numeric,
            default = 10
        ),
        # Additional SC-specific metadata
        cell_metadata = new_computed_property(
            function(self) {
                if (!requireNamespace("Seurat", quietly = TRUE)) {
                    stop("Package 'Seurat' must be installed to use this function.")
                }
                meta <- self@seurat@meta.data
                meta$cell_id <- rownames(meta)
                # Filter cells that have the required identity
                meta <- meta[!is.na(meta[[self@cell_identity]]), ]
                return(meta)
            },
            name = "cell_metadata"
        )
    )
)

#' @title Convert Seurat Object to Single-Cell PhyloExpressionSet
#' @description Convert a Seurat object with phylostratum information into a 
#' ScPhyloExpressionSet object for single-cell phylotranscriptomic analysis.
#' 
#' @param seurat A Seurat object containing single-cell expression data
#' @param phylomap A data frame with two columns: phylostratum assignments and gene IDs
#' @param cell_identity Character string specifying which cell identity to use from Seurat metadata
#' @param slot Character string specifying which slot to use from the Seurat object (default: "data")
#' @param name A character string naming the dataset (default: "Single-Cell Phylo Expression Set")
#' @param min_cells_per_type Minimum number of cells required per cell type (default: 10)
#' @param strata_labels Optional character vector of labels for phylostrata
#' @param ... Additional arguments passed to ScPhyloExpressionSet constructor
#' 
#' @return A ScPhyloExpressionSet object
#' 
#' @examples
#' # Convert Seurat object to ScPhyloExpressionSet
#' # sc_phyex_set <- as_ScPhyloExpressionSet(seurat_obj, phylo_map, "celltype",
#' #                                        name = "Brain Development SC Dataset")
#' 
#' @importFrom dplyr inner_join
#' @export
as_ScPhyloExpressionSet <- function(seurat, 
                                   phylomap,
                                   cell_identity,
                                   slot = "data",
                                   name = "Single-Cell Phylo Expression Set",
                                   min_cells_per_type = 10,
                                   strata_labels = NULL,
                                   ...) {
    
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }
    
    # Validate inputs
    if (!cell_identity %in% colnames(seurat@meta.data)) {
        stop(paste("Cell identity", cell_identity, "not found in Seurat metadata"))
    }
    
    # Get raw single-cell expression data
    sc_counts <- .get_sc_expression_matrix(seurat, slot)
    
    # Match with phylomap
    colnames(phylomap) <- c("Stratum", "GeneID")
    gene_ids <- rownames(sc_counts)
    
    # Filter phylomap to genes present in data
    phylomap_filtered <- phylomap[phylomap$GeneID %in% gene_ids, ]
    
    # Filter expression data to genes with phylostratum info
    common_genes <- intersect(gene_ids, phylomap_filtered$GeneID)
    sc_counts <- sc_counts[common_genes, , drop = FALSE]
    phylomap_filtered <- phylomap_filtered[phylomap_filtered$GeneID %in% common_genes, ]
    
    # Order phylomap to match gene order in expression data
    phylomap_ordered <- phylomap_filtered[match(rownames(sc_counts), phylomap_filtered$GeneID), ]
    
    # Create strata factor
    if (is.null(strata_labels)) {
        strata_labels <- sort(unique(as.numeric(phylomap_ordered$Stratum)))
    }
    strata <- factor(as.numeric(phylomap_ordered$Stratum), 
                    levels = sort(unique(as.numeric(phylomap_ordered$Stratum))), 
                    labels = strata_labels)
    names(strata) <- rownames(sc_counts)
    
    # Get cell identities and filter
    cell_types <- seurat@meta.data[[cell_identity]]
    names(cell_types) <- colnames(seurat)
    
    # Filter to cells that have valid identities and are in our expression matrix
    valid_cells <- !is.na(cell_types) & names(cell_types) %in% colnames(sc_counts)
    sc_counts <- sc_counts[, names(cell_types)[valid_cells], drop = FALSE]
    cell_types <- cell_types[valid_cells]
    
    # Filter cell types with sufficient cells for pseudobulking
    cell_type_counts <- table(cell_types)
    valid_cell_types <- names(cell_type_counts)[cell_type_counts >= min_cells_per_type]
    
    if (length(valid_cell_types) == 0) {
        stop(paste("No cell types have at least", min_cells_per_type, "cells"))
    }
    
    # Keep only cells from valid cell types
    final_valid_cells <- cell_types %in% valid_cell_types
    sc_counts <- sc_counts[, final_valid_cells, drop = FALSE]
    cell_types <- cell_types[final_valid_cells]
    
    # Create groups factor for each cell
    groups <- factor(cell_types, levels = valid_cell_types)
    
    return(ScPhyloExpressionSet(
        seurat = seurat,
        cell_identity = cell_identity,
        slot = slot,
        min_cells_per_type = min_cells_per_type,
        strata = strata,
        gene_ids = rownames(sc_counts),
        counts = sc_counts,  # Raw single-cell data
        groups = groups,     # Cell-level groups
        name = name,
        index_type = "TAI",
        conditions_label = "Cell Types",
        is_time_series = FALSE,
        ...
    ))
}

#' @title Get Single-Cell Expression Matrix
#' @description Internal function to extract expression matrix from Seurat object.
#' 
#' @param seurat A Seurat object
#' @param slot Character string specifying which slot to use
#' 
#' @return Expression matrix with genes as rows and cells as columns
#' 
#' @keywords internal
.get_sc_expression_matrix <- function(seurat, slot) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }
    
    # Get expression data
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
    
    return(as.matrix(counts))
}

#' @title Pseudobulk Single-Cell Expression Data
#' @description Internal function to create pseudobulk expression data from single-cell data
#' by aggregating cells within each cell type/identity.
#' 
#' @param seurat A Seurat object
#' @param cell_identity Character string specifying which cell identity to use
#' @param slot Character string specifying which slot to use
#' @param min_cells_per_type Minimum number of cells required per cell type
#' 
#' @return Matrix with pseudobulked expression data, one column per cell type
#' 
#' @keywords internal
.pseudobulk_expression <- function(seurat, cell_identity, slot, min_cells_per_type) {
    # Get raw single-cell data
    sc_counts <- .get_sc_expression_matrix(seurat, slot)
    
    # Get cell identities
    cell_types <- seurat@meta.data[[cell_identity]]
    names(cell_types) <- colnames(seurat)
    
    # Remove cells with NA identity
    valid_cells <- !is.na(cell_types) & names(cell_types) %in% colnames(sc_counts)
    sc_counts <- sc_counts[, names(cell_types)[valid_cells], drop = FALSE]
    cell_types <- cell_types[valid_cells]
    
    # Filter cell types with sufficient cells
    cell_type_counts <- table(cell_types)
    valid_cell_types <- names(cell_type_counts)[cell_type_counts >= min_cells_per_type]
    
    if (length(valid_cell_types) == 0) {
        stop(paste("No cell types have at least", min_cells_per_type, "cells"))
    }
    
    # Keep only cells from valid cell types
    valid_cells_final <- cell_types %in% valid_cell_types
    sc_counts <- sc_counts[, valid_cells_final, drop = FALSE]
    cell_types <- cell_types[valid_cells_final]
    
    # Pseudobulk by summing expression within each cell type
    pseudobulk_list <- lapply(valid_cell_types, function(ct) {
        ct_cells <- names(cell_types)[cell_types == ct]
        if (length(ct_cells) == 1) {
            as.matrix(sc_counts[, ct_cells, drop = FALSE])
        } else {
            Matrix::rowSums(sc_counts[, ct_cells, drop = FALSE])
        }
    })
    
    pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
    colnames(pseudobulk_matrix) <- valid_cell_types
    
    return(as.matrix(pseudobulk_matrix))
}

#' @title Print Method for ScPhyloExpressionSet
#' @description Print summary information for a ScPhyloExpressionSet object.
S7::method(print, ScPhyloExpressionSet) <- function(x, ...) {
    cat("ScPhyloExpressionSet object\n")
    cat("Name:", x@name, "\n")
    cat("Species:", ifelse(is.null(x@species), "Not specified", x@species), "\n")
    cat("Index type:", x@index_type, "\n")
    cat("Cell identity used:", x@cell_identity, "\n")
    cat("Expression slot used:", x@slot, "\n")
    cat("Min cells per type:", x@min_cells_per_type, "\n")
    cat("\nDimensions:\n")
    cat("- Genes:", x@num_genes, "\n")
    cat("- Cell types:", x@num_conditions, "\n")
    cat("- Phylostrata:", x@num_strata, "\n")
    cat("- Total cells:", ncol(x@seurat), "\n")
    cat("\nCell types:\n")
    print(table(x@cell_metadata[[x@cell_identity]]))
}

# S7 method extensions for ScPhyloExpressionSet

#' @export
S7::method(collapse, ScPhyloExpressionSet) <- function(phyex_set) {
    # For single-cell data, collapse returns pseudobulked data
    # So we create a regular PhyloExpressionSet from the pseudobulked data
    phylomap <- data.frame(
        Stratum = as.numeric(phyex_set@strata),
        GeneID = phyex_set@gene_ids
    )
    
    pseudobulk_data <- cbind(
        phylomap,
        phyex_set@counts_collapsed  # Use inherited counts_collapsed property
    )
    
    as_PhyloExpressionSet(
        data = pseudobulk_data,
        groups = colnames(phyex_set@counts_collapsed),
        name = paste(phyex_set@name, "collapsed"),
        species = phyex_set@species,
        index_type = phyex_set@index_type,
        conditions_label = phyex_set@conditions_label,
        is_time_series = phyex_set@is_time_series
    )
}

