#' @title Single-Cell PhyloExpressionSet Class
#' @description S7 class for single-cell phylotranscriptomic expression data.
#' This class handles Seurat objects and provides pseudobulking functionality.
#' 
#' @param strata Factor vector of phylostratum assignments for each gene
#' @param strata_values Numeric vector of phylostratum values used in TXI calculations
#' @param gene_ids Character vector of gene identifiers
#' @param name Character string naming the dataset (default: "Phylo Expression Set")
#' @param species Character string specifying the species (default: NULL)
#' @param index_type Character string specifying the transcriptomic index type (default: "TXI")
#' @param identities_label Character string labeling the identities (default: "Cell Type")
#' @param null_conservation_sample_size Numeric value for null conservation sample size (default: 5000)
#' @param precomputed_null_conservation_txis Precomputed null conservation TXI values (default: NULL)
#' @param seurat A Seurat object containing single-cell expression data
#' @param layer Character string specifying which layer to use from the Seurat object (default: "data")
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
        layer = new_property(
            class = class_character,
            default = "counts"
        ),
        
        identities_label = new_property(
            class = class_character,
            default = "Cell Type"
        ),
        
        ## IMPLEMENTED ABSTRACT PROPERTIES
        expression = new_property(
            getter = function(self) .get_expression_matrix(self@seurat, self@layer)
        ),
        
        expression_collapsed = new_property(
            getter = function(self) .pseudobulk_expression(self@seurat, self@layer)
        ),
        groups = new_property(
            class = class_factor,
            getter = function(self) Seurat::Idents(self@seurat)
        ),
        
        ## SC-SPECIFIC TXI
        TXI_sample = new_property(
            class = class_double,
            getter = function(self) .TXI_sc(self@expression, self@strata_values)
        ),
        cell_metadata = new_property(
            getter = function(self) self@seurat@meta.data
        )
    )
)

## CONSTRUCTORS

#' @title Convert Seurat Object to Single-Cell PhyloExpressionSet
#' @description Convert a Seurat object with phylostratum information into a 
#' ScPhyloExpressionSet object for single-cell phylotranscriptomic analysis.
#' 
#' @param seurat A Seurat object containing single-cell expression data
#' @param strata Factor vector of phylostratum assignments for each gene
#' @param layer Character string specifying which layer to use from the Seurat object (default: "counts")
#' @param name A character string naming the dataset (default: "Single-Cell Phylo Expression Set")
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
                                    strata,
                                    layer = "counts",
                                    name = "Single-Cell Phylo Expression Set",
                                    ...) {
    
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }

    # Get raw single-cell expression data
    sc_counts <- .get_expression_matrix(seurat, layer)

    names(strata) <- rownames(sc_counts)
    strata_values <- as.numeric(strata)
    names(strata_values) <- rownames(sc_counts)
    
    
    
    return(ScPhyloExpressionSet(
        seurat = seurat,
        layer = layer,
        strata = strata,
        strata_values = strata_values,
        gene_ids = rownames(sc_counts),
        name = name,
        ...
    ))
}

#' @title Match Single-Cell Expression Data with Phylostratum Map
#' @description Join single-cell gene expression data (from a Seurat object) with a phylostratum mapping to create 
#' a ScPhyloExpressionSet object.
#' 
#' @param seurat A Seurat object containing single-cell expression data
#' @param phylomap A data frame with two columns: phylostratum assignments and gene IDs
#' @param layer Character string specifying which layer to use from the Seurat object (default: "counts")
#' @param strata_legend A data frame with two columns: phylostratum assignments and name of each stratum. If NULL, no labels will be added (default: NULL)
#' If NULL, uses sorted unique values from column 1
#' @param ... Additional arguments passed to as_ScPhyloExpressionSet
#' 
#' @return A ScPhyloExpressionSet object
#' 
#' @examples
#' # Match Seurat object with phylostratum map
#' # sc_set <- match_map_sc(seurat_obj, phylo_map, layer = "counts", name = "SC Matched Dataset")
#' 
#' @importFrom dplyr inner_join relocate
#' @export
match_map_sc <- function(seurat, 
                         phylomap,
                         layer = "counts",
                         strata_legend = NULL,
                         ...) {

    # Match with phylomap
    colnames(phylomap) <- c("Stratum", "GeneID")

    sc_counts <- .get_expression_matrix(seurat, layer)
    gene_ids <- rownames(sc_counts)

    # Filter phylomap to genes present in data and remove duplicates
    phylomap_filtered <- phylomap[phylomap$GeneID %in% gene_ids, ]
    phylomap_filtered <- phylomap_filtered[!duplicated(phylomap_filtered$GeneID), ]

    # Filter expression data to genes with phylostratum info
    common_genes <- intersect(gene_ids, phylomap_filtered$GeneID)

    # Subset sc_counts and Seurat object to common genes
    sc_counts <- sc_counts[common_genes, , drop = FALSE]
    seurat <- subset(seurat, features = common_genes)

    # Order phylomap to match gene order in expression data
    matched_idx <- match(rownames(sc_counts), phylomap_filtered$GeneID)
    phylomap_ordered <- phylomap_filtered[matched_idx, ]

    # Now create strata
    
    if (is.null(strata_legend)) {
        levels <- sort(unique(as.numeric(phylomap$Stratum)))
        labels <- levels
    }
    else {
        levels <- strata_legend[[1]]
        labels <- strata_legend[[2]]
    }
    strata <- factor(as.numeric(phylomap_ordered$Stratum), 
                     levels = levels, 
                     labels = labels)

    names(strata) <- rownames(sc_counts)

    return(as_ScPhyloExpressionSet(seurat=seurat, strata = strata, layer = layer, ...))
}

## HELPER FUNCTIONS

#' @title Get Expression Matrix from Seurat Object
#' @description Extract expression matrix from Seurat object, preserving sparse format when possible.
#' 
#' @param seurat A Seurat object
#' @param layer Character string specifying which layer to use
#' @return Expression matrix (sparse or dense)
#' 
#' @keywords internal
.get_expression_matrix <- function(seurat, layer) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("Package 'SeuratObject' must be installed to use this function.")
    }
    expr <- SeuratObject::LayerData(seurat, layer = layer)
    # Always return a 2D sparse matrix, even for single gene/cell
    if (is.null(dim(expr)) || length(dim(expr)) == 1) {
        if (!requireNamespace("Matrix", quietly = TRUE)) {
            stop("Package 'Matrix' must be installed to use this function.")
        }
        # Try to infer if this is a single gene (row) or single cell (column)
        if (!is.null(names(expr))) {
            # If names match rownames, it's a single cell; if names match colnames, it's a single gene
            if (!is.null(rownames(seurat)) && all(names(expr) %in% rownames(seurat))) {
                # Single cell, multiple genes
                expr <- Matrix::Matrix(expr, nrow = length(expr), ncol = 1, sparse = TRUE,
                                       dimnames = list(names(expr), colnames(seurat)))
            } else if (!is.null(colnames(seurat)) && all(names(expr) %in% colnames(seurat))) {
                # Single gene, multiple cells
                expr <- Matrix::Matrix(expr, nrow = 1, ncol = length(expr), sparse = TRUE,
                                       dimnames = list(rownames(seurat), names(expr)))
            } else {
                # Fallback: treat as column vector
                expr <- Matrix::Matrix(expr, nrow = length(expr), ncol = 1, sparse = TRUE)
            }
        } else {
            # Single value
            expr <- Matrix::Matrix(expr, nrow = 1, ncol = 1, sparse = TRUE,
                                   dimnames = list(rownames(seurat), colnames(seurat)))
        }
    }
    return(expr)
}

#' @title Create Pseudobulk Expression Data
#' @description Aggregate single-cell data by cell identity (using Seurat Idents) to create pseudobulk expression.
#' 
#' @param seurat A Seurat object
#' @param layer Expression layer to use
#' @return Matrix of pseudobulked expression
#' 
#' @keywords internal
.pseudobulk_expression <- function(seurat, layer) {
    # Get expression data
    counts <- .get_expression_matrix(seurat, layer)
    
    # Get cell identities from Seurat Idents (factor: names = cell names, values = cell types)
    cell_types <- Seurat::Idents(seurat)
    
    # Use all unique identities present (levels of the factor)
    unique_cell_types <- levels(cell_types)
    
    # Pseudobulk by summing within each cell type
    pseudobulk_list <- lapply(unique_cell_types, function(ct) {
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
    
    # If only one cell per type, ensure output is a matrix
    result <- do.call(cbind, pseudobulk_list)
    colnames(result) <- unique_cell_types
    
    # Convert to regular matrix for consistency with base class
    return(as.matrix(result))
}

#' @title Calculate TXI for single cell expression sparse matrix.
#' @description Internal function to calculate TXI for expression data.
#' 
#' @param expression Matrix of expression values, dgmatrix
#' @param strata_values Numeric vector of phylostratum values
#' @return Vector of TXI values
#' 
#' @keywords internal
.TXI_sc <- function(expression, strata_values) {
    # Calculate column sums
    col_sums <- Matrix::colSums(expression)
    
    # Handle zero column sums (cells with no expression)
    zero_cols <- col_sums == 0
    
    if (all(zero_cols)) {
        # If all cells have zero expression, return vector of NAs
        txi <- rep(NA_real_, ncol(expression))
    } else {
        # Calculate TXI only for non-zero columns
        txi <- rep(NA_real_, ncol(expression))
        
        if (any(!zero_cols)) {
            non_zero_expr <- expression[, !zero_cols, drop = FALSE]
            non_zero_sums <- col_sums[!zero_cols]
            
            txi_non_zero <- as.numeric((Matrix::t(non_zero_expr) %*% strata_values) / non_zero_sums)
            txi[!zero_cols] <- txi_non_zero
        }
    }
    
    names(txi) <- colnames(expression)
    return(txi)
}

## METHOD IMPLEMENTATIONS

#' @export
S7::method(collapse, ScPhyloExpressionSet) <- function(phyex_set) {
    # For single-cell data, collapse returns pseudobulked data
    # So we create a regular BulkPhyloExpressionSet from the pseudobulked data
    phylomap <- data.frame(
        Stratum = phyex_set@strata_values,
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

    if (length(valid_indices) < length(gene_indices))
        warning("Some of the specified genes were not found in the dataset")

    if (length(valid_indices) == 0) {
        stop("None of the specified genes were found in the dataset")
    }

    if (length(valid_indices) < 2) {
        stop("Please provide at least 2 genes for select_genes() on a ScPhyloExpressionSet. Single-gene selection is not supported.")
    }

    genes <- phyex_set@gene_ids[valid_indices]

    seurat_subset <- subset(phyex_set@seurat, features = genes)
    as_ScPhyloExpressionSet(
        seurat = seurat_subset,
        strata = phyex_set@strata[valid_indices], 
        layer = phyex_set@layer,
        name = phyex_set@name,
        species = phyex_set@species,
        index_type = phyex_set@index_type
    )
}

## PRINT METHODS

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
    cat("Expression layer used:", x@layer, "\n")
    cat("Total cells:", ncol(x@seurat), "\n")
    cat("Valid cells:", x@num_samples, "\n")
    cat("Cells per type:\n")
    print(table(x@groups))
}


#' @title Downsample Expression Matrix
#' @description Downsample the cells in a ScPhyloExpressionSet and return the expression matrix as a dense matrix.
#' @param phyex_set A ScPhyloExpressionSet object
#' @param downsample Integer, number of cells to keep per identity (default: 10)
#' @return A dense expression matrix (genes x downsampled cells)
#' @details
#' This function randomly downsamples the cells in the Seurat object of a ScPhyloExpressionSet
#' and returns the resulting expression matrix as a regular R matrix (not sparse).
#' Useful for quick plotting or prototyping with large single-cell datasets.
#' @examples
#' # Downsample to 20 cells per identity and get the matrix
#' # mat <- downsample_expression(sc_phyex_set, downsample = 20)
#' @export
downsample_expression <- function(phyex_set, downsample = 10) {
    seurat <- subset(phyex_set@seurat, downsample = downsample)
    as.matrix(.get_expression_matrix(seurat, phyex_set@layer))
}

#' @title Get available identities
#' @description Return the available metadata columns that can be used as identities in the Seurat object.
#' @param phyex_set A ScPhyloExpressionSet object
#' @return Character vector of available identity columns
#' @export
available_identities <- function(phyex_set) {
    colnames(phyex_set@seurat@meta.data)
}

#' @title Set Identities for ScPhyloExpressionSet
#' @description Change the Seurat identities to a different metadata column.
#' @param phyex_set A ScPhyloExpressionSet object
#' @param identity_name Character, name of the metadata column to use as new identities
#' @return ScPhyloExpressionSet object with updated identities
#' @export
set_identities <- function(phyex_set, identity_name) {
    available <- available_identities(phyex_set)
    if (!(identity_name %in% available)) {
        stop(
            sprintf(
                "Column '%s' not found in Seurat metadata. Available options are: %s",
                identity_name, paste(available, collapse = ", ")
            )
        )
    }
    Idents(phyex_set@seurat) <- identity_name
    phyex_set@identities_label <- identity_name
    phyex_set
}

#' @title Set Identity Order for ScPhyloExpressionSet
#' @description Change the order of the identities (factor levels) in the Seurat object.
#' @param phyex_set A ScPhyloExpressionSet object
#' @param new_order Character vector specifying the new order of identities
#' @return ScPhyloExpressionSet object with updated identity order
#' @export
reorder_identities <- function(phyex_set, new_order) {
    current_idents <- Seurat::Idents(phyex_set@seurat)
    if (!all(new_order %in% levels(current_idents))) {
        stop(
            sprintf(
                "Some values in new_order are not current identities. Current identities are: %s",
                paste(levels(current_idents), collapse = ", ")
            )
        )
    }
    Seurat::Idents(phyex_set@seurat) <- factor(current_idents, levels = new_order)
    phyex_set
}

#' @title Check if object is a ScPhyloExpressionSet
#' @description Checks if the input is a PhyloExpressionSet S7 object and throws an error if not.
#' @param phyex_set An object to check
#' @return Invisibly returns TRUE if check passes, otherwise throws an error
#' @export
check_ScPhyloExpressionSet <- function(phyex_set) {
    if (!S7::S7_inherits(phyex_set, ScPhyloExpressionSet)) {
        stop("Input must be a ScPhyloExpressionSet S7 object.", call. = FALSE)
    }
    invisible(TRUE)
}