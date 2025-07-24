#' @title Single-Cell PhyloExpressionSet Class
#' @description S7 class for single-cell phylotranscriptomic expression data.
#' This class handles Seurat objects and provides pseudobulking functionality.
#' 
#' @param strata Factor vector of phylostratum assignments for each gene
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
            getter = function(self) .TXI_sc(self@expression, self@strata)
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
#' @param phylomap A data frame with two columns: phylostratum assignments and gene IDs
#' @param layer Character string specifying which layer to use from the Seurat object (default: "counts")
#' @param name A character string naming the dataset (default: "Single-Cell Phylo Expression Set")
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
                                    layer = "counts",
                                    name = "Single-Cell Phylo Expression Set",
                                    strata_labels = NULL,
                                    ...) {
    
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }

    # Get raw single-cell expression data
    sc_counts <- .get_expression_matrix(seurat, layer)
    
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
        layer = layer,
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
#' @param expression_matrix Matrix of expression values, dgmatrix
#' @param strata Vector of phylostratum assignments
#' @return Vector of TXI values
#' 
#' @keywords internal
.TXI_sc <- function(expression, strata) {
    txi <- (Matrix::t(expression) %*% strata) / Matrix::colSums(expression) |>
        as.vector()
    names(txi) <- colnames(expression)
    return(txi)
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
            expr_matrix <- .get_expression_matrix(phyex_set@seurat, phyex_set@layer)
            expr_subset <- expr_matrix[selected_genes_in_seurat, , drop = FALSE]
            
            # Create new Seurat object with subset data
            Seurat::CreateSeuratObject(
                counts = if (phyex_set@layer == "counts") expr_subset else NULL,
                data = if (phyex_set@layer == "data") expr_subset else NULL,
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
        layer = phyex_set@layer,
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
    cat("Expression layer used:", x@layer, "\n")
    cat("Total cells:", ncol(x@seurat), "\n")
    cat("Valid cells:", x@num_samples, "\n")
    cat("Cells per type:\n")
    print(table(x@groups))
}
