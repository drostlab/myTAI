#' @title Single-Cell PhyloExpressionSet Class
#' @description S7 class for single-cell phylotranscriptomic expression data.
#' This class stores expression matrices and metadata, with support for dimensional
#' reductions and pseudobulking functionality.
#' 
#' @param strata Factor vector of phylostratum assignments for each gene
#' @param strata_values Numeric vector of phylostratum values used in TXI calculations
#' @param expression Sparse or dense matrix of expression counts with genes as rows and cells as columns
#' @param groups Factor vector indicating which identity each cell belongs to (derived from selected_idents column in metadata)
#' @param name Character string naming the dataset (default: "Phylo Expression Set")
#' @param species Character string specifying the species (default: NULL)
#' @param index_type Character string specifying the transcriptomic index type (default: "TXI")
#' @param identities_label Character string labeling the identities (default: "Cell Type")
#' @param metadata Data frame with cell metadata, where rownames correspond to cell IDs and columns contain cell attributes
#' @param selected_idents Character string specifying which metadata column is currently used for grouping cells
#' @param reductions List of dimensional reduction matrices (PCA, UMAP, etc.) with cells as rows and dimensions as columns
#' @param idents_colours List of named character vectors specifying colors for each identity level, organized by metadata column name
#' @param null_conservation_sample_size Numeric value for null conservation sample size (default: 5000)
#' @param .null_conservation_txis Precomputed null conservation TXI values (default: NULL)
#' @param .pseudobulk_cache Internal cache for pseudobulked expression matrices by different groupings
#' @param .TXI_sample Internal storage for computed TXI values
#' 
#' @details
#' The ScPhyloExpressionSet class provides a comprehensive framework for single-cell 
#' phylotranscriptomic analysis. Key features include:
#' 
#' \strong{Identity Management:}
#' The \code{selected_idents} property determines which metadata column is used for grouping cells.
#' When changed, it automatically updates the \code{groups} property and invalidates cached
#' pseudobulk data to ensure consistency.
#' 
#' \strong{Dimensional Reductions:}
#' The \code{reductions} property stores pre-computed dimensional reductions (PCA, UMAP, etc.).
#' If not provided during construction from Seurat objects, basic PCA and UMAP are computed
#' automatically.
#' 
#' \strong{Color Management:}
#' \code{idents_colours} allows custom color schemes for different metadata columns, ensuring
#' consistent visualization across plots.
#' 
#' \strong{Computed Properties:}
#' Several properties are computed automatically when accessed:
#' \itemize{
#'   \item \code{available_idents} - Character vector of factor columns in metadata that can be used for grouping (automatically detected from metadata)
#'   \item \code{expression_collapsed} - Matrix of pseudobulked expression data (genes x identities), created by summing expression within each identity group
#'   \item \code{TXI_sample} - Named numeric vector of TXI (Transcriptomic Age Index) values for each cell, computed using efficient C++ implementation
#' }
#' 
#' Inherited computed properties from PhyloExpressionSetBase include:
#' \itemize{
#'   \item \code{gene_ids} - Character vector of gene identifiers
#'   \item \code{identities} - Character vector of identity labels  
#'   \item \code{sample_names} - Character vector of sample names (cell IDs)
#'   \item \code{num_identities} - Integer count of unique cell types/identities
#'   \item \code{num_samples} - Integer count of total cells
#'   \item \code{num_genes} - Integer count of genes
#'   \item \code{num_strata} - Integer count of phylostrata
#'   \item \code{index_full_name} - Full name of the transcriptomic index type
#'   \item \code{group_map} - List mapping identity names to cell IDs
#'   \item \code{TXI} - Numeric vector of TXI values for each identity (computed from pseudobulked expression)
#'   \item \code{null_conservation_txis} - Matrix of null conservation TXI values for statistical testing
#' }
#' These properties use lazy evaluation and caching for optimal performance.
#' 
#' @examples
#' \dontrun{
#' # Create from Seurat object
#' sc_set <- ScPhyloExpressionSet_from_seurat(seurat_obj, strata)
#' 
#' # Switch to different cell grouping
#' sc_set@selected_idents <- "development_stage"
#' 
#' # Access pseudobulked data (computed automatically)
#' pseudobulk <- sc_set@expression_collapsed
#' 
#' # Access TXI values for each cell
#' txi_values <- sc_set@TXI_sample
#' }
#' 
#' @import S7
#' @export
ScPhyloExpressionSet <- new_class("ScPhyloExpressionSet",
    parent = PhyloExpressionSetBase,
    properties = list(
        expression = new_required_property(
            # class = class_matrix (S7 doesn't have a class matrix yet)
            validator = function(value) {
                if (any(is.na(value))) return("cannot contain NA values. Check expression matrix.")
                if (length(value) == 0) return("cannot be empty. Check expression matrix.")
            },
            setter = function(self, value) {
                if (!length(value))
                    return(self)
                # recompute TXI sample
                self@.TXI_sample <- .TXI_sc_adaptive(value, self@strata_values)
                self@.pseudobulk_cache <- list()
                self@expression <- value
                self
            },
            name = "expression"
        ),
        .pseudobulk_cache = new_property(
            class = class_list,
            default = list()
        ),
        expression_collapsed = new_property(
            getter = function(self) {
                key <- self@selected_idents
                cache <- self@.pseudobulk_cache
                if (!is.null(cache[[key]]))
                    return(cache[[key]])
                pb <- .pseudobulk_expression(self@expression, self@groups)
                rownames(pb) <- self@gene_ids
                colnames(pb) <- levels(self@groups)
                cache[[key]] <- pb
                self@.pseudobulk_cache <- cache
                return(pb)
            }
        ),
        .TXI_sample = new_property(
            class = class_double,
            default = NULL
        ),
        TXI_sample = new_property(
            class = class_double,
            validator = function(value) {
                if (any(is.na(value))) return("cannot contain NA values. Check expression matrix.")
                if (length(value) == 0) return("cannot be empty. Check expression matrix.")
            },
            getter = \(self) self@.TXI_sample
        ),
        identities_label = new_property(
            class = class_character,
            default = "Cell Type"
        ),
        metadata = new_property(
            class = class_data.frame,
            
            setter = function(self, value) {
                if (!length(value))
                    return(self)
                
                # Convert character and discrete numeric columns to factors for consistent handling
                for (col_name in colnames(value)) {
                    col <- value[[col_name]]
                    if (!is.factor(col)) {
                        if (is.character(col) || is.integer(col)) {
                            # Always convert character columns to factors
                            unique_vals <- unique(col)
                            value[[col_name]] <- factor(col, levels = unique_vals)
                        } 
                    }
                }
                self@metadata <- value
                self
            },
            default = NULL
        ),
        available_idents = new_property(
            class = class_character,
            getter = function(self) {
                if (is.null(self@metadata))
                    return(character())
                cols <- colnames(self@metadata)
                discrete_cols <- cols[sapply(self@metadata, \(col) is.factor(col))]
                discrete_cols
            }
        ),
        selected_idents = new_property(
            class = class_character,
            setter = function(self, value) {
                if (!length(value)) 
                    return(self)
                # update the groups property
                if (!(value %in% self@available_idents))
                    stop("@selected_idents must be a factor column in self@metadata")
                self@groups <- self@metadata[[value]]
                self@selected_idents <- value
                self
            }
        ),
        idents_colours = new_property(
            class = class_list,
            default = list(),
            validator = function(value) {
                if (is.null(value) || length(value) == 0) return()
                
                # Check that each element is a named character vector
                for (i in seq_along(value)) {
                    element <- value[[i]]
                    if (!is.character(element)) {
                        return(sprintf(" element %d must be a character vector", i))
                    }
                    if (is.null(names(element)) || any(names(element) == "")) {
                        return(sprintf(" element %d must be a named character vector with all names non-empty", i))
                    }
                }
            }
        ),
        # dimensionality reductions
        reductions = new_property(
            class = class_list,
            default = list(),
            validator = function(value) {
                if (is.null(value) || length(value) == 0) return()
                
                # Check that each element is a matrix with proper dimensions
                for (reduction_name in names(value)) {
                    reduction <- value[[reduction_name]]
                    if (!is.matrix(reduction) && !is.data.frame(reduction)) {
                        return(sprintf("Reduction '%s' must be a matrix or data.frame", reduction_name))
                    }
                    if (ncol(reduction) < 2) {
                        return(sprintf("Reduction '%s' must have at least 2 dimensions", reduction_name))
                    }
                }
            }
        )
    ),
    validator = function(self) {
        # Validator for idents_colours property
        if (!is.null(self@idents_colours) && length(self@idents_colours) > 0) {
            # Check that names correspond to available identities
            invalid_names <- setdiff(names(self@idents_colours), self@available_idents)
            if (length(invalid_names) > 0) {
                return(
                    sprintf(
                        "@idents_colours must have its names be a subset of available_idents: %s. Invalid: %s",
                        paste(self@available_idents, collapse = ", "),
                        paste(invalid_names, collapse = ", ")
                    )
                )
            }
            
            # Check that all identity values are covered and nothing more
            for (identity_name in names(self@idents_colours)) {
                if (is.null(self@metadata) || !(identity_name %in% colnames(self@metadata))) {
                    return(sprintf("Identity '%s' not found in metadata", identity_name))
                }
                
                colours <- self@idents_colours[[identity_name]]
                identity_values <- unique(as.character(self@metadata[[identity_name]]))
                colour_names <- names(colours)
                
                # Check for missing identity values
                missing_values <- setdiff(identity_values, colour_names)
                if (length(missing_values) > 0) {
                    return(
                        sprintf(
                            "idents_colours for '%s' missing colours for values: %s",
                            identity_name,
                            paste(missing_values, collapse = ", ")
                        )
                    )
                }
                
                # Check for extra colour names not in identity values
                extra_colours <- setdiff(colour_names, identity_values)
                if (length(extra_colours) > 0) {
                    return(
                        sprintf(
                            "idents_colours for '%s' has colours for non-existent values: %s",
                            identity_name,
                            paste(extra_colours, collapse = ", ")
                        )
                    )
                }
            }
        }
        
        # Validator for .pseudobulk_cache property
        if (!is.null(self@.pseudobulk_cache) && length(self@.pseudobulk_cache) > 0) {
            cache_keys <- names(self@.pseudobulk_cache)
            if (!is.null(cache_keys) && length(cache_keys) > 0) {
                # Check that cache keys are a subset of available identities
                invalid_cache_keys <- setdiff(cache_keys, self@available_idents)
                if (length(invalid_cache_keys) > 0) {
                    return(
                        sprintf(
                            "@.pseudobulk_cache keys must be a subset of available_idents: %s. Invalid: %s",
                            paste(self@available_idents, collapse = ", "),
                            paste(invalid_cache_keys, collapse = ", ")
                        )
                    )
                }
            }
        }

        # Validate metadata
        if (!is.null(self@metadata))
            if (nrow(self@metadata) != self@num_samples)
                return("@metadata must have the same number of rows as the number of samples")

        # Validate expression_collapsed rownames match gene_ids
        # if (!identical(rownames(self@expression_collapsed), self@gene_ids)) {
        #     return("@expression_collapsed rownames must match @gene_ids")
        # }
        
        # unique_groups <- sort(as.character(unique(self@groups)))
        # if (!identical(sort(self@identities), unique_groups)) {
        #     return("@identities must match unique values in @groups")
        # }
    }
)

## CONSTRUCTORS

#' @title Convert Seurat Object to Single-Cell PhyloExpressionSet
#' @description Convert a Seurat object with phylostratum information into a 
#' ScPhyloExpressionSet object for single-cell phylotranscriptomic analysis.
#' Automatically extracts dimensional reductions if present, or computes basic
#' PCA and UMAP if none are available.
#' 
#' @param seurat A Seurat object containing single-cell expression data
#' @param strata Factor vector of phylostratum assignments for each gene
#' @param layer Character string specifying which layer to use from the Seurat object (default: "counts")
#' @param selected_idents Character string specifying which metadata column to use for grouping (default: NULL, uses active idents)
#' @param name A character string naming the dataset (default: "Single-Cell Phylo Expression Set")
#' @param seed Integer seed for reproducible UMAP computation (default: 42)
#' @param ... Additional arguments passed to ScPhyloExpressionSet constructor
#' 
#' @return A ScPhyloExpressionSet object
#' 
#' @examples
#' # Convert Seurat object to ScPhyloExpressionSet
#' # sc_set <- ScPhyloExpressionSet_from_seurat(seurat_obj, strata,
#' #                                           name = "Brain Development SC Dataset")
#' 
#' @export
ScPhyloExpressionSet_from_seurat <- function(seurat, 
                                             strata,
                                             layer = "counts",
                                             selected_idents = NULL,
                                             name = "Single-Cell PhyloExpressionSet",
                                             seed = 42,
                                             ...) {
    
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }

    # Get raw single-cell expression data
    sc_counts <- .get_expression_from_seurat(seurat, layer)
    
    # Extract metadata from Seurat object
    metadata <- seurat@meta.data
    
    # Get current identities and add to metadata if not already present
    if (is.null(selected_idents)) {
        groups <- Seurat::Idents(seurat)
        metadata$active_ident <- groups
        selected_idents <- "active_ident"
    } else {
        if (!(selected_idents %in% colnames(metadata))) {
            stop(sprintf("selected_idents '%s' must be a column name in seurat@meta.data", selected_idents))
        }
        groups <- metadata[[selected_idents]]
    }

    names(strata) <- rownames(sc_counts)
    strata_values <- as.numeric(strata)
    names(strata_values) <- rownames(sc_counts)
    
    # Extract dimensional reductions from Seurat object
    reductions <- list()
    if (length(seurat@reductions) > 0) {
        for (reduction_name in names(seurat@reductions)) {
            embeddings <- Seurat::Embeddings(seurat, reduction = reduction_name)
            reductions[[reduction_name]] <- embeddings
        }
    }
    
    # If no reductions are available, compute basic ones
    if (length(reductions) == 0) {
        message("No dimensional reductions found in Seurat object. Computing PCA and UMAP...")
        
        # Compute PCA
        pca_coords <- .compute_reduction(sc_counts, method = "PCA", seed = seed)
        reductions[["pca"]] <- pca_coords
        
        # Compute UMAP (if uwot is available)
        if (requireNamespace("uwot", quietly = TRUE)) {
            umap_coords <- .compute_reduction(sc_counts, method = "UMAP", seed = seed)
            reductions[["umap"]] <- umap_coords
        }
    }
    
    return(ScPhyloExpressionSet(
        strata = strata,
        strata_values = strata_values,
        expression = sc_counts,
        groups = groups,
        name = name,
        metadata = metadata,
        .TXI_sample = .TXI_sc_adaptive(sc_counts, strata_values),
        selected_idents = selected_idents,
        reductions = reductions,
        ...
    ))
}

#' @title Create Single-Cell PhyloExpressionSet from Expression Matrix
#' @description Create a ScPhyloExpressionSet object from an expression matrix and metadata.
#' 
#' @param expression_matrix Sparse or dense expression matrix with genes as rows and cells as columns
#' @param strata Factor vector of phylostratum assignments for each gene
#' @param metadata Data frame with cell metadata, rownames should match colnames of expression_matrix
#' @param groups_column Character string specifying which metadata column to use for initial grouping (default: first factor column found)
#' @param name A character string naming the dataset (default: "Single-Cell Phylo Expression Set")
#' @param ... Additional arguments passed to ScPhyloExpressionSet constructor
#' 
#' @details
#' This function creates a ScPhyloExpressionSet from basic components. The \code{groups_column}
#' parameter determines the initial \code{selected_idents} value, which can be changed later
#' using the setter. All discrete columns in metadata are automatically converted to factors
#' for consistent handling.
#' 
#' @return A ScPhyloExpressionSet object
#' 
#' @examples
#' \dontrun{
#' # Create from matrix and metadata
#' sc_set <- ScPhyloExpressionSet_from_matrix(expr_matrix, strata, metadata,
#'                                           groups_column = "cell_type")
#' 
#' # Later change grouping
#' sc_set@selected_idents <- "treatment_group"
#' }
#' 
#' @export
ScPhyloExpressionSet_from_matrix <- function(expression_matrix,
                                             strata,
                                             metadata,
                                             groups_column = NULL,
                                             name = "Single-Cell Phylo Expression Set",
                                             ...) {
    
    # Validate inputs
    if (nrow(expression_matrix) != length(strata)) {
        stop("Number of genes in expression matrix must match length of strata vector")
    }
    
    if (ncol(expression_matrix) != nrow(metadata)) {
        stop("Number of cells in expression matrix must match number of rows in metadata")
    }
    
    # Set names if missing
    names(strata) <- rownames(expression_matrix)
    strata_values <- as.numeric(strata)
    names(strata_values) <- rownames(expression_matrix)
    
    # Determine groups column
    if (is.null(groups_column)) {
        # Find first factor column in metadata
        factor_cols <- sapply(metadata, is.factor)
        if (any(factor_cols)) {
            groups_column <- names(factor_cols)[which(factor_cols)[1]]
            message("Using '", groups_column, "' as groups column")
        } else {
            stop("No factor columns found in metadata. Please specify groups_column or ensure metadata has factor columns.")
        }
    }
    
    if (!(groups_column %in% colnames(metadata))) {
        stop("groups_column '", groups_column, "' not found in metadata")
    }
    
    groups <- metadata[[groups_column]]
    selected_idents <- groups_column

    
    return(ScPhyloExpressionSet(
        strata = strata,
        strata_values = strata_values,
        expression = expression_matrix,
        groups = groups,
        name = name,
        metadata = metadata,
        .TXI_sample = .TXI_sc_adaptive(expression_matrix, strata_values),
        selected_idents = selected_idents,
        ...
    ))
}

#' @title Match Single-Cell Expression Data with Phylostratum Map (Seurat)
#' @description Join single-cell gene expression data (from a Seurat object) with a phylostratum mapping to create 
#' a ScPhyloExpressionSet object. Automatically extracts dimensional reductions and metadata.
#' 
#' @param seurat A Seurat object containing single-cell expression data
#' @param phylomap A data frame with two columns: phylostratum assignments and gene IDs
#' @param layer Character string specifying which layer to use from the Seurat object (default: "counts")
#' @param strata_legend A data frame with two columns: phylostratum assignments and name of each stratum. If NULL, numeric labels will be used (default: NULL)
#' @param selected_idents Character string specifying which metadata column to use for initial grouping (default: NULL, uses active idents)
#' @param seed Integer seed for reproducible UMAP computation when reductions need to be computed (default: 42)
#' @param ... Additional arguments passed to ScPhyloExpressionSet_from_seurat
#' 
#' @details
#' This is a convenience function that combines phylostratum mapping with Seurat object conversion.
#' Only genes present in both the expression data and phylomap will be retained. The function
#' extracts all metadata and dimensional reductions from the Seurat object.
#' 
#' @return A ScPhyloExpressionSet object
#' 
#' @examples
#' \dontrun{
#' # Match Seurat object with phylostratum map
#' sc_set <- match_map_sc_seurat(seurat_obj, phylo_map, layer = "counts")
#' 
#' # With custom grouping column
#' sc_set <- match_map_sc_seurat(seurat_obj, phylo_map, 
#'                              selected_idents = "development_stage")
#' }
#' 
#' @export
match_map_sc_seurat <- function(seurat, 
                                phylomap,
                                layer = "counts",
                                strata_legend = NULL,
                                ...) {

    # Match with phylomap
    colnames(phylomap) <- c("Stratum", "GeneID")

    sc_counts <- .get_expression_from_seurat(seurat, layer)
    gene_ids <- rownames(sc_counts)

    # Filter phylomap to genes present in data and remove duplicates
    phylomap_filtered <- phylomap[phylomap$GeneID %in% gene_ids, ]
    phylomap_filtered <- phylomap_filtered[!duplicated(phylomap_filtered$GeneID), ]

    # Filter expression data to genes with phylostratum info
    common_genes <- intersect(gene_ids, phylomap_filtered$GeneID)

    # Subset sc_counts to common genes
    sc_counts <- sc_counts[common_genes, , drop = FALSE]
    
    # Subset seurat object to common genes for extracting reductions
    seurat <- subset(seurat, features = common_genes)

    # Order phylomap to match gene order in expression data
    matched_idx <- match(rownames(sc_counts), phylomap_filtered$GeneID)
    phylomap_ordered <- phylomap_filtered[matched_idx, ]

    # Create strata
    if (is.null(strata_legend)) {
        levels <- sort(unique(as.numeric(phylomap$Stratum)))
        labels <- levels
    } else {
        levels <- strata_legend[[1]]
        labels <- strata_legend[[2]]
    }
    strata <- factor(as.numeric(phylomap_ordered$Stratum), 
                     levels = levels, 
                     labels = labels)

    names(strata) <- rownames(sc_counts)

    return(ScPhyloExpressionSet_from_seurat(seurat = seurat, strata = strata, layer = layer, ...))
}

#' @title Match Expression Matrix with Phylostratum Map
#' @description Join single-cell gene expression matrix with a phylostratum mapping to create 
#' a ScPhyloExpressionSet object.
#' 
#' @param expression_matrix Expression matrix with genes as rows and cells as columns
#' @param metadata Data frame with cell metadata, rownames should match colnames of expression_matrix
#' @param phylomap A data frame with two columns: phylostratum assignments and gene IDs
#' @param strata_legend A data frame with two columns: phylostratum assignments and name of each stratum. If NULL, numeric labels will be used (default: NULL)
#' @param groups_column Character string specifying which metadata column to use for initial grouping (default: first factor column found)
#' @param name Character string naming the dataset (default: derived from function)
#' @param ... Additional arguments passed to ScPhyloExpressionSet_from_matrix
#' 
#' @details
#' This function combines phylostratum mapping with expression matrix and metadata to create
#' a ScPhyloExpressionSet. Only genes present in both the expression matrix and phylomap
#' will be retained. All discrete metadata columns are converted to factors automatically.
#' 
#' @return A ScPhyloExpressionSet object
#' 
#' @examples
#' \dontrun{
#' # Match expression matrix with phylostratum map
#' sc_set <- match_map_sc_matrix(expr_matrix, metadata, phylo_map)
#' 
#' # With specific grouping column
#' sc_set <- match_map_sc_matrix(expr_matrix, metadata, phylo_map,
#'                              groups_column = "cell_type")
#' }
#' 
#' @export
match_map_sc_matrix <- function(expression_matrix,
                                metadata, 
                                phylomap,
                                strata_legend = NULL,
                                ...) {

    # Match with phylomap
    colnames(phylomap) <- c("Stratum", "GeneID")

    gene_ids <- rownames(expression_matrix)

    # Filter phylomap to genes present in data and remove duplicates
    phylomap_filtered <- phylomap[phylomap$GeneID %in% gene_ids, ]
    phylomap_filtered <- phylomap_filtered[!duplicated(phylomap_filtered$GeneID), ]

    # Filter expression data to genes with phylostratum info
    common_genes <- intersect(gene_ids, phylomap_filtered$GeneID)

    # Subset expression matrix to common genes
    expression_matrix <- expression_matrix[common_genes, , drop = FALSE]

    # Order phylomap to match gene order in expression data
    matched_idx <- match(rownames(expression_matrix), phylomap_filtered$GeneID)
    phylomap_ordered <- phylomap_filtered[matched_idx, ]

    # Create strata
    if (is.null(strata_legend)) {
        levels <- sort(unique(as.numeric(phylomap$Stratum)))
        labels <- levels
    } else {
        levels <- strata_legend[[1]]
        labels <- strata_legend[[2]]
    }
    strata <- factor(as.numeric(phylomap_ordered$Stratum), 
                     levels = levels, 
                     labels = labels)

    names(strata) <- rownames(expression_matrix)

    return(ScPhyloExpressionSet_from_matrix(
        expression_matrix = expression_matrix, 
        strata = strata, 
        metadata = metadata,
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
.get_expression_from_seurat <- function(seurat, layer) {
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
#' @description Aggregate single-cell data by cell identity to create pseudobulk expression.
#' 
#' @param expression Expression matrix with genes as rows and cells as columns
#' @param groups Factor vector indicating which group each cell belongs to
#' @return Matrix of pseudobulked expression
#' 
#' @keywords internal
.pseudobulk_expression <- function(expression, groups) {
    # Use all unique identities present (levels of the factor)
    unique_groups <- levels(groups)
    
    # Pseudobulk by summing within each group
    pseudobulk_list <- lapply(unique_groups, function(group) {
        group_cells <- which(groups == group)
        if (length(group_cells) == 1) {
            expression[, group_cells, drop = FALSE]
        } else {
            if (inherits(expression, "sparseMatrix")) {
                Matrix::rowSums(expression[, group_cells, drop = FALSE])
            } else {
                rowSums(expression[, group_cells, drop = FALSE])
            }
        }
    })
    
    # If only one cell per type, ensure output is a matrix
    result <- do.call(cbind, pseudobulk_list)
    colnames(result) <- unique_groups
    
    # Convert to regular matrix for consistency with base class
    return(as.matrix(result))
}

#' @title Compute Dimensional Reduction
#' @description Compute PCA or UMAP on expression data when not available in stored reductions.
#' 
#' @param expression Expression matrix with genes as rows and cells as columns
#' @param method Character string: "PCA" or "UMAP"
#' @param seed Integer seed for reproducible results (default: 42)
#' @return Matrix with cells as rows and dimensions as columns
#' 
#' @keywords internal
.compute_reduction <- function(expression, method = c("PCA", "UMAP"), seed = 42) {
    method <- match.arg(method)
    
    # Prepare expression data
    expr <- log1p(as.matrix(expression))
    # Remove genes with zero variance
    nonzero_var_genes <- apply(expr, 1, function(x) stats::var(x) != 0)
    expr <- expr[nonzero_var_genes, , drop = FALSE]
    
    if (method == "PCA") {
        coords <- stats::prcomp(t(expr), scale. = TRUE)$x[, 1:2]
    } else if (method == "UMAP") {
        if (!requireNamespace("uwot", quietly = TRUE)) {
            stop("Package 'uwot' must be installed to compute UMAP.")
        }
        set.seed(seed)
        coords <- uwot::umap(t(expr), scale = TRUE)
        rownames(coords) <- colnames(expression)
    }
    
    return(coords)
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

#' @title Adaptive TXI calculation for single cell expression
#' @description Automatically selects the best TXI implementation based on dataset size
#' and available computational resources. Uses R implementation for smaller datasets
#' and C++ implementation with optimal parallelization for larger datasets.
#' 
#' @param expression Matrix of expression values, dgCMatrix
#' @param strata_values Numeric vector of phylostratum values
#' @param force_method Character, force specific method: "r", "cpp_simple", or "cpp_batched"
#' @param ncores Integer, number of cores to use (default: parallel::detectCores())
#' @return Vector of TXI values
#' 
#' @details Based on performance benchmarking:
#' - R implementation is fastest for < 50,000 cells
#' - C++ batched implementation becomes advantageous for >= 100,000 cells
#' - Optimal core count scales with dataset size
#' 
#' @keywords internal
.TXI_sc_adaptive <- function(expression, strata_values, force_method = NULL, ncores = NULL) {
    n_cells <- ncol(expression)
    n_genes <- nrow(expression)
    
    # Detect available cores if not specified
    if (is.null(ncores)) {
        ncores <- min(parallel::detectCores(), 16)
    }
    
    # Determine optimal method based on benchmarking results
    if (!is.null(force_method)) {
        method <- force_method
    } else if (n_cells < 100000) {
        # R implementation dominates for smaller datasets
        method <- "r"
    } else {
        # C++ batched with high parallelization for very large datasets
        method <- "cpp_batched"
        # Optimal core count for large datasets
        ncores <- min(ncores, max(8, ncores))
    }
    
    # Execute with selected method
    switch(method,
           "r" = .TXI_sc(expression, strata_values),
           "cpp_batched" = {
               # Adaptive batch size based on dataset size
               batch_size <- max(2000, min(5000, n_cells %/% ncores))
               cpp_txi_sc(expression, strata_values, batch_size = batch_size, ncores = ncores)
           },
           stop("Unknown method: ", method)
    )
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
        groups = phyex_set@groups,
        name = paste(phyex_set@name, "collapsed"),
        species = phyex_set@species,
        index_type = phyex_set@index_type,
        identities_label = phyex_set@identities_label
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

    obj <- S7::valid_eventually(phyex_set, function(x) {
        x@strata <- phyex_set@strata[valid_indices]
        x@strata_values <- phyex_set@strata_values[valid_indices]
        x@expression <- phyex_set@expression[valid_indices, , drop = FALSE]
        x
    })
    
    obj
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
    cat("Total cells:", x@num_samples, "\n")
    cat("Cells per type:\n")
    print(table(x@groups))
    
    # Print metadata information
    if (!is.null(x@metadata) && ncol(x@metadata) > 0) {
        cat("Available metadata:\n")
        for (col_name in x@available_idents) {
            col_data <- x@metadata[[col_name]]
            unique_vals <- levels(col_data)
            cat("  ", col_name, ": ", paste(unique_vals, collapse = ", "), "\n", sep = "")
        }
    }
}


#' @title Downsample Expression Matrix by Groups
#' @description Downsample an expression matrix by randomly selecting a specified number of samples from each group.
#' @param expression_matrix Expression matrix with genes as rows and samples as columns
#' @param groups Factor vector indicating which group each sample belongs to
#' @param downsample Integer, number of samples to keep per group (default: 10)
#' @return A dense expression matrix (genes x downsampled samples)
#' @details
#' This function randomly samples up to \code{downsample} samples from each group level.
#' The returned expression matrix is converted to dense format and maintains column names
#' from the original matrix. Useful for creating balanced subsets for visualization or
#' when memory is limited.
#' @examples
#' \dontrun{
#' # Downsample expression matrix to 5 samples per group
#' downsampled <- downsample_expression(expr_matrix, groups, downsample = 5)
#' }
#' @export
downsample_expression <- function(expression_matrix, groups, downsample = 10) {
    # Get unique groups
    unique_groups <- levels(groups)
    
    # Sample cells from each group
    sampled_indices <- c()
    for (group in unique_groups) {
        group_indices <- which(groups == group)
        n_sample <- min(downsample, length(group_indices))
        if (n_sample > 0) {
            sampled_indices <- c(sampled_indices, sample(group_indices, n_sample))
        }
    }
    
    # Return downsampled expression matrix as dense matrix
    result <- as.matrix(expression_matrix[, sampled_indices, drop = FALSE])
    return(result)
}

#' @title Downsample ScPhyloExpressionSet
#' @description Create a downsampled copy of a ScPhyloExpressionSet object with fewer cells per identity.
#' @param phyex_set A ScPhyloExpressionSet object
#' @param downsample Integer, number of cells to keep per identity (default: 10)
#' @return A new ScPhyloExpressionSet object with downsampled cells
#' @details
#' This function creates a new ScPhyloExpressionSet with a subset of cells, maintaining
#' the same proportional representation across identities. The sampling is stratified 
#' by the current \code{selected_idents} grouping. All metadata and reductions are 
#' filtered to match the selected cells.
#' @examples
#' \dontrun{
#' # Downsample to 20 cells per identity
#' small_set <- downsample(sc_phyex_set, downsample = 20)
#' 
#' # Change grouping and downsample
#' sc_phyex_set@selected_idents <- "treatment"
#' treatment_set <- downsample(sc_phyex_set, downsample = 15)
#' }
#' @export
downsample <- function(phyex_set, downsample = 10) {
    # Use downsample_expression to get the downsampled matrix
    downsampled_expr <- downsample_expression(phyex_set@expression, phyex_set@groups, downsample)
    
    # Get the indices of selected cells
    selected_cell_names <- colnames(downsampled_expr)
    selected_indices <- match(selected_cell_names, colnames(phyex_set@expression))
    
    # Filter metadata to selected cells
    filtered_metadata <- phyex_set@metadata[selected_indices, , drop = FALSE]
    
    # Filter groups to selected cells
    filtered_groups <- phyex_set@groups[selected_indices]
    
    # Filter reductions to selected cells if they exist
    filtered_reductions <- list()
    if (length(phyex_set@reductions) > 0) {
        for (reduction_name in names(phyex_set@reductions)) {
            reduction <- phyex_set@reductions[[reduction_name]]
            # Match by rownames (cell names)
            reduction_cell_indices <- match(selected_cell_names, rownames(reduction))
            valid_indices <- !is.na(reduction_cell_indices)
            if (any(valid_indices)) {
                filtered_reductions[[reduction_name]] <- reduction[reduction_cell_indices[valid_indices], , drop = FALSE]
            }
        }
    }
    
    # Create new ScPhyloExpressionSet with downsampled data
    new_phyex_set <- ScPhyloExpressionSet(
        strata = phyex_set@strata,
        strata_values = phyex_set@strata_values,
        expression = downsampled_expr,
        groups = filtered_groups,
        name = paste(phyex_set@name, "downsampled"),
        species = phyex_set@species,
        index_type = phyex_set@index_type,
        identities_label = phyex_set@identities_label,
        metadata = filtered_metadata,
        selected_idents = phyex_set@selected_idents,
        reductions = filtered_reductions,
        idents_colours = phyex_set@idents_colours,
        null_conservation_sample_size = phyex_set@null_conservation_sample_size
    )
    
    # Warn about subsampling
    total_cells <- length(phyex_set@groups)
    n_selected <- ncol(downsampled_expr)
    if (n_selected < total_cells) {
        message(paste("Downsampled from", total_cells, "to", n_selected, "cells."))
    }
    
    return(new_phyex_set)
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