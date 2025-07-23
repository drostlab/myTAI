#' @title Single-Cell PhyloExpressionSet Class
#' @description S7 class for single-cell phylotranscriptomic expression data.
#' This class handles Seurat objects and provides pseudobulking functionality.
#' 
#' @param strata Factor vector of phylostratum assignments for each gene
#' @param gene_ids Character vector of gene identifiers
#' @param name Character string naming the dataset (default: "Phylo Expression Set")
#' @param species Character string specifying the species (default: NULL)
#' @param index_type Character string specifying the transcriptomic index type (default: "TXI")
#' @param identities_label Character string labeling the identities (default: "Cell Types")
#' @param null_conservation_sample_size Numeric value for null conservation sample size (default: 5000)
#' @param precomputed_null_conservation_txis Precomputed null conservation TXI values (default: NULL)
#' @param seurat A Seurat object containing single-cell expression data
#' @param cell_identity Character string specifying which cell identity to use from Seurat metadata
#' @param slot Character string specifying which slot to use from the Seurat object (default: "data")
#' @param min_cells_per_identity Minimum number of cells required per cell type (default: 10)
#' 
#' @import S7
#' @export
ScPhyloExpressionSet <- new_class("ScPhyloExpressionSet",
    parent = PhyloExpressionSetBase,
    properties = list(
        ## SC-SPECIFIC REQUIRED PROPERTIES
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
        min_cells_per_identity = new_property(
            class = class_numeric,
            default = 10
        ),
        
        ## RAW EXPRESSION DATA
        expression = new_property(
            getter = function(self) .get_expression_matrix(self@seurat, self@slot)
        ),
        
        ## IMPLEMENTED ABSTRACT PROPERTIES
        identities = new_property(
            class = class_factor,
            getter = function(self) .get_valid_cell_types(self@seurat, self@cell_identity, self@min_cells_per_identity)
        ),
        identities_label = new_property(
            class = class_character,
            default = "Cell Types"
        ),
        expression_collapsed = new_property(
            getter = function(self) .pseudobulk_expression(self@seurat, self@cell_identity, self@slot, self@min_cells_per_identity)
        ),
        
        ## SC-SPECIFIC TXI
        TXI_sample = new_property(
            class = class_double,
            getter = function(self) {
                # Efficient per-cell TXI calculation
                counts <- .get_expression_matrix(self@seurat, self@slot)
                strata_numeric <- as.numeric(self@strata)
                txi <- (Matrix::t(counts) %*% strata_numeric) / Matrix::colSums(counts)
                txi_vector <- as.vector(txi)
                names(txi_vector) <- colnames(counts)
                return(txi_vector)
            }
        ),
        
        ## SC-SPECIFIC PROPERTIES
        sample_names = new_property(
            class = class_character,
            getter = function(self) {
                # Get cell identities and filter to valid ones
                cell_types <- self@seurat@meta.data[[self@cell_identity]]
                names(cell_types) <- colnames(self@seurat)
                valid_cells <- cell_types %in% self@identities
                return(names(cell_types)[valid_cells])
            }
        ),
        groups = new_property(
            class = class_factor,
            getter = function(self) .map_cells_to_groups(self@seurat, self@cell_identity, self@identities)
        ),
        cell_metadata = new_property(
            getter = function(self) .get_cell_metadata(self@seurat, self@cell_identity)
        )
    )
)

## CONSTRUCTORS

#' @title Convert Seurat Object to Single-Cell PhyloExpressionSet
#' @description Convert a Seurat object with phylostratum information into a 
#' ScPhyloExpressionSet object for single-cell phylotranscriptomic analysis.
#' 
#' @param seurat A Seurat object containing single-cell expression data
#' @param phylomap A data frame with two columns: phylostratum assignments and gene IDs
#' @param cell_identity Character string specifying which cell identity to use from Seurat metadata
#' @param slot Character string specifying which slot to use from the Seurat object (default: "data")
#' @param name A character string naming the dataset (default: "Single-Cell Phylo Expression Set")
#' @param min_cells_per_identity Minimum number of cells required per cell type (default: 10)
#' @param strata_labels Optional character vector of labels for phylostrata
#' @param ... Additional arguments passed to ScPhyloExpressionSet constructor
#' 
#' @return A ScPhyloExpressionSet object
#' 
#' @examples
#' # Convert Seurat object to ScPhyloExpressionSet
#' # sc_set <- as_ScPhyloExpressionSet(seurat_obj, phylo_map, "celltype",
#' #                                  name = "Brain Development SC Dataset")
#' 
#' @importFrom dplyr inner_join
#' @export
as_ScPhyloExpressionSet <- function(seurat, 
                                   phylomap,
                                   cell_identity,
                                   slot = "data",
                                   name = "Single-Cell Phylo Expression Set",
                                   min_cells_per_identity = 10,
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
    sc_counts <- .get_expression_matrix(seurat, slot)
    
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
    
    return(ScPhyloExpressionSet(
        seurat = seurat,
        cell_identity = cell_identity,
        slot = slot,
        min_cells_per_identity = min_cells_per_identity,
        strata = strata,
        gene_ids = rownames(sc_counts),
        name = name,
        ...
    ))
}

## HELPER FUNCTIONS

#' @title Get Expression Matrix from Seurat Object
#' @description Extract expression matrix from Seurat object, preserving sparse format when possible.
#' 
#' @param seurat A Seurat object
#' @param slot Character string specifying which slot to use
#' @return Expression matrix (sparse or dense)
#' 
#' @keywords internal
.get_expression_matrix <- function(seurat, slot) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }
    
    counts <- tryCatch({
        # Try new API first (layer parameter)
        if (slot == "data") {
            # Try to get data layer, fall back to counts if empty
            tryCatch({
                result <- Seurat::GetAssayData(seurat, layer = "data")
                if (nrow(result) == 0 || ncol(result) == 0) {
                    Seurat::GetAssayData(seurat, layer = "counts")
                } else {
                    result
                }
            }, error = function(e) {
                Seurat::GetAssayData(seurat, layer = "counts")
            })
        } else if (slot == "counts") {
            Seurat::GetAssayData(seurat, layer = "counts")
        } else if (slot == "scale.data") {
            tryCatch({
                Seurat::GetAssayData(seurat, layer = "scale.data")
            }, error = function(e) {
                # Fall back to old API for scale.data
                suppressWarnings(Seurat::GetAssayData(seurat, slot = "scale.data"))
            })
        } else {
            # For any other slot name, try as layer first, then as slot
            tryCatch({
                Seurat::GetAssayData(seurat, layer = slot)
            }, error = function(e) {
                suppressWarnings(Seurat::GetAssayData(seurat, slot = slot))
            })
        }
    }, error = function(e) {
        # Final fallback to old API
        suppressWarnings(Seurat::GetAssayData(seurat, slot = slot))
    })
    
    return(counts)
}

#' @title Get Valid Cell Types
#' @description Get cell types that meet minimum cell count requirements.
#' 
#' @param seurat A Seurat object
#' @param cell_identity Column name in metadata
#' @param min_cells Minimum number of cells per type
#' @return Factor of valid cell types
#' 
#' @keywords internal
.get_valid_cell_types <- function(seurat, cell_identity, min_cells) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }
    
    cell_types <- seurat@meta.data[[cell_identity]]
    cell_type_counts <- table(cell_types)
    valid_types <- names(cell_type_counts)[cell_type_counts >= min_cells]
    
    if (length(valid_types) == 0) {
        stop(paste("No cell types have at least", min_cells, "cells"))
    }
    
    return(factor(valid_types, levels = valid_types, ordered = TRUE))
}

#' @title Create Pseudobulk Expression Data
#' @description Aggregate single-cell data by cell type to create pseudobulk expression.
#' 
#' @param seurat A Seurat object
#' @param cell_identity Column name in metadata
#' @param slot Expression slot to use
#' @param min_cells Minimum cells per type
#' @return Matrix of pseudobulked expression
#' 
#' @keywords internal
.pseudobulk_expression <- function(seurat, cell_identity, slot, min_cells) {
    # Get expression data
    counts <- .get_expression_matrix(seurat, slot)
    
    # Get cell identities
    cell_types <- seurat@meta.data[[cell_identity]]
    names(cell_types) <- colnames(seurat)
    
    # Filter to valid cells and types
    valid_cells <- !is.na(cell_types) & names(cell_types) %in% colnames(counts)
    counts <- counts[, names(cell_types)[valid_cells], drop = FALSE]
    cell_types <- cell_types[valid_cells]
    
    # Get valid cell types
    cell_type_counts <- table(cell_types)
    valid_cell_types <- names(cell_type_counts)[cell_type_counts >= min_cells]
    
    if (length(valid_cell_types) == 0) {
        stop(paste("No cell types have at least", min_cells, "cells"))
    }
    
    # Keep only cells from valid types
    valid_cells_final <- cell_types %in% valid_cell_types
    counts <- counts[, valid_cells_final, drop = FALSE]
    cell_types <- cell_types[valid_cells_final]
    
    # Pseudobulk by summing within each cell type
    pseudobulk_list <- lapply(valid_cell_types, function(ct) {
        ct_cells <- names(cell_types)[cell_types == ct]
        if (length(ct_cells) == 1) {
            counts[, ct_cells, drop = FALSE]
        } else {
            if (inherits(counts, "sparseMatrix")) {
                Matrix::rowSums(counts[, ct_cells, drop = FALSE])
            } else {
                rowSums(counts[, ct_cells, drop = FALSE])
            }
        }
    })
    
    result <- do.call(cbind, pseudobulk_list)
    colnames(result) <- valid_cell_types
    
    # Convert to regular matrix for consistency with base class
    return(as.matrix(result))
}

#' @title Map Cells to Groups
#' @description Create factor mapping each cell to its group/identity.
#' 
#' @param seurat A Seurat object
#' @param cell_identity Column name in metadata
#' @param valid_identities Valid identity levels
#' @return Factor mapping cells to groups
#' 
#' @keywords internal
.map_cells_to_groups <- function(seurat, cell_identity, valid_identities) {
    cell_types <- seurat@meta.data[[cell_identity]]
    names(cell_types) <- colnames(seurat)
    
    # Filter to valid identities and maintain order
    valid_cells <- cell_types %in% valid_identities
    cell_types_filtered <- cell_types[valid_cells]
    
    # Return factor with same order as filtered cell names
    return(factor(cell_types_filtered, levels = valid_identities))
}

#' @title Get Cell Metadata
#' @description Extract relevant cell metadata from Seurat object.
#' 
#' @param seurat A Seurat object
#' @param cell_identity Column name for cell identity
#' @return Data frame of cell metadata
#' 
#' @keywords internal
.get_cell_metadata <- function(seurat, cell_identity) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }
    
    meta <- seurat@meta.data
    meta$cell_id <- rownames(meta)
    
    # Filter cells that have the required identity
    meta <- meta[!is.na(meta[[cell_identity]]), ]
    
    return(meta)
}

## METHOD IMPLEMENTATIONS

#' @export
S7::method(collapse, ScPhyloExpressionSet) <- function(phyex_set) {
    # For single-cell data, collapse returns pseudobulked data
    # So we create a regular BulkPhyloExpressionSet from the pseudobulked data
    phylomap <- data.frame(
        Stratum = as.numeric(phyex_set@strata),
        GeneID = phyex_set@gene_ids
    )
    
    pseudobulk_data <- cbind(
        phylomap,
        phyex_set@expression_collapsed  # Use inherited expression_collapsed property
    )
    
    as_BulkPhyloExpressionSet(
        data = pseudobulk_data,
        groups = colnames(phyex_set@expression_collapsed),
        name = paste(phyex_set@name, "collapsed"),
        species = phyex_set@species,
        index_type = phyex_set@index_type
    )
}

#' @export
S7::method(select_genes, ScPhyloExpressionSet) <- function(phyex_set, genes) {
    # Find indices of selected genes
    gene_indices <- match(genes, phyex_set@gene_ids)
    
    # Remove any NA indices (genes not found)
    valid_indices <- gene_indices[!is.na(gene_indices)]
    
    if (length(valid_indices) == 0) {
        stop("None of the specified genes were found in the dataset")
    }
    
    # Create new seurat object with selected genes
    selected_genes_in_seurat <- intersect(genes, rownames(phyex_set@seurat))
    
    if (length(selected_genes_in_seurat) == 0) {
        stop("None of the specified genes were found in the Seurat object")
    }
    
    # Subset the Seurat object using Seurat's subset function
    seurat_subset <- tryCatch({
        # Try using Seurat's subset function
        Seurat::subset(phyex_set@seurat, features = selected_genes_in_seurat)
    }, error = function(e) {
        # Fallback: create new Seurat object with subset data
        tryCatch({
            # Extract the expression matrix for selected genes
            expr_matrix <- .get_expression_matrix(phyex_set@seurat, phyex_set@slot)
            expr_subset <- expr_matrix[selected_genes_in_seurat, , drop = FALSE]
            
            # Create new Seurat object with subset data
            Seurat::CreateSeuratObject(
                counts = if (phyex_set@slot == "counts") expr_subset else NULL,
                data = if (phyex_set@slot == "data") expr_subset else NULL,
                meta.data = phyex_set@seurat@meta.data
            )
        }, error = function(e2) {
            # Final fallback: just use the original seurat object
            # This is not ideal but allows the function to continue
            warning("Could not subset Seurat object, using original object")
            phyex_set@seurat
        })
    })
    
    # Create phylomap for selected genes
    phylomap <- data.frame(
        Stratum = as.numeric(phyex_set@strata[valid_indices]),
        GeneID = phyex_set@gene_ids[valid_indices]
    )
    
    as_ScPhyloExpressionSet(
        seurat = seurat_subset,
        phylomap = phylomap,
        cell_identity = phyex_set@cell_identity,
        slot = phyex_set@slot,
        min_cells_per_identity = phyex_set@min_cells_per_identity,
        name = phyex_set@name,
        species = phyex_set@species,
        index_type = phyex_set@index_type
    )
}

## PRINT METHODS

#' @export
S7::method(print, ScPhyloExpressionSet) <- function(x, ...) {
    # Print base information (inline parent method)
    cat("PhyloExpressionSet object\n")
    cat("Class:", class(x)[[1]], "\n")
    cat("Name:", x@name, "\n")
    cat("Species:", ifelse(is.null(x@species), "Not specified", x@species), "\n")
    cat("Index type:", x@index_type, "\n")
    cat(x@identities_label, ":", paste(as.character(x@identities), collapse = ", "), "\n")
    cat("Number of genes:", x@num_genes, "\n")
    cat("Number of", tolower(x@identities_label), ":", x@num_identities, "\n")
    cat("Number of phylostrata:", x@num_strata, "\n")
    
    # Print single-cell specific information
    cat("Cell identity used:", x@cell_identity, "\n")
    cat("Expression slot used:", x@slot, "\n")
    cat("Min cells per type:", x@min_cells_per_identity, "\n")
    cat("Total cells:", ncol(x@seurat), "\n")
    cat("Cells per type:\n")
    print(table(x@cell_metadata[[x@cell_identity]]))
}
