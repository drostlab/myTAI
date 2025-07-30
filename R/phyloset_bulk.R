#' @title Bulk PhyloExpressionSet Class
#' @description S7 class for bulk RNA-seq phylotranscriptomic expression data.
#' This class handles expression data with biological replicates.
#' 
#' @param strata Factor vector of phylostratum assignments for each gene
#' @param gene_ids Character vector of gene identifiers
#' @param .expression Matrix of expression counts with genes as rows and samples as columns
#' @param .groups Factor vector indicating which identity each sample belongs to
#' @param name Character string naming the dataset (default: "Phylo Expression Set")
#' @param species Character string specifying the species (default: NULL)
#' @param index_type Character string specifying the transcriptomic index type (default: "TXI")
#' @param identities_label Character string labeling the identities (default: "Identities")
#' @param null_conservation_sample_size Numeric value for null conservation sample size (default: 5000)
#' @param precomputed_null_conservation_txis Precomputed null conservation TXI values (default: NULL)
#' 
#' @import S7
#' @export
BulkPhyloExpressionSet <- new_class("BulkPhyloExpressionSet",
    parent = PhyloExpressionSetBase,
    properties = list(
        ## BULK-SPECIFIC REQUIRED PROPERTIES
        .expression = new_required_property(
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check expression data."
                if (length(value) == 0) "cannot be empty. Check expression data."
            },
            name = ".expression"
        ),
        .groups = new_required_property(
            class = class_factor,
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check groups."
                if (length(value) == 0) "cannot be empty. Check groups."
            },
            name = ".groups"
        ),
        
        ## IMPLEMENTED ABSTRACT PROPERTIES
        expression = new_property(
            getter = function(self) self@.expression
        ),
        expression_collapsed = new_property(
            getter = function(self) .collapse_replicates(self@.expression, self@.groups)
        ),
        groups = new_property(
            getter = function(self) self@.groups
        )
    )
)


#' @title Convert Data to BulkPhyloExpressionSet
#' @description Convert a data frame with phylostratum, gene ID, and expression data 
#' into a BulkPhyloExpressionSet object.
#' 
#' @param data A data frame where column 1 contains phylostratum information, 
#' column 2 contains gene IDs, and columns 3+ contain expression data
#' @param groups A factor or character vector indicating which group each sample belongs to.
#' Default uses column names from expression data
#' @param name A character string naming the dataset. Default uses the variable name
#' @param strata_legend A data frame with two columns: phylostratum assignments and name of each stratum. If NULL, no labels will be added (default: NULL)
#' If NULL, uses sorted unique values from column 1
#' @param ... Additional arguments passed to BulkPhyloExpressionSet constructor
#' 
#' @return A BulkPhyloExpressionSet object
#' 
#' @examples
#' # Convert data frame to BulkPhyloExpressionSet
#' # bulk_set <- as_BulkPhyloExpressionSet(my_data, 
#' #                                      groups = c("stage1", "stage1", "stage2", "stage2"),
#' #                                      name = "Development Dataset")
#' 
#' @export
as_BulkPhyloExpressionSet <- function(data, 
                                      groups = colnames(data[, 3:ncol(data)]),
                                      name = deparse(substitute(data)),
                                      strata_legend = NULL,
                                      ...) {
    gene_ids <- as.character(data[[2]])

    if (is.null(strata_legend)) {
        levels <- sort(unique(as.numeric(data[[1]])))
        labels <- levels
    }
    else {
        levels <- strata_legend[[1]]
        labels <- strata_legend[[2]]
    }
    strata <- factor(as.numeric(data[[1]]), levels=levels, labels=labels)

    names(strata) <- gene_ids

    groups <- factor(groups, levels=unique(groups))
    
    expression <- as.matrix(data[3:ncol(data)])
    rownames(expression) <- gene_ids
    
    obj <- BulkPhyloExpressionSet(
        strata = strata,
        gene_ids = gene_ids,
        .expression = expression,
        .groups = groups,
        name = name,
        ...
    )


    return(obj)
}

#' @title Match Gene Expression Data with Phylostratum Map
#' @description Join gene expression data with a phylostratum mapping to create 
#' a BulkPhyloExpressionSet object.
#' 
#' @param data A data frame where column 1 contains gene IDs and columns 2+ contain expression data
#' @param phylomap A data frame with two columns: phylostratum assignments and gene IDs
#' @param groups A factor or character vector indicating which group each sample belongs to.
#' Default uses column names from expression data
#' @param name A character string naming the dataset. Default uses the variable name
#' @param ... Additional arguments passed to as_BulkPhyloExpressionSet
#' 
#' @return A BulkPhyloExpressionSet object
#' 
#' @examples
#' # Match expression data with phylostratum map
#' # bulk_set <- match_map(expression_data, phylo_map, 
#' #                       groups = c("stage1", "stage2", "stage3"),
#' #                       name = "Matched Dataset")
#' 
#' @importFrom dplyr inner_join relocate
#' @export
match_map <- function(data, 
                      phylomap,
                      groups = colnames(data[,2:ncol(data)]),
                      name = NULL,
                      ...) {
    if (is.null(name)) name <- deparse(substitute(data))
    colnames(phylomap) <- c("Stratum", "GeneID")
    
    data <- data |>
        inner_join(phylomap, by="GeneID") |>
        relocate(Stratum, .before = 1)
    
    return(as_BulkPhyloExpressionSet(data, groups = groups, name = name, ...))
}


#' @title Collapse Expression Data Across Replicates
#' @description Internal function to collapse expression data across replicates by taking row means.
#' 
#' @param expression Matrix of expression counts
#' @param groups Factor indicating group membership
#' @return Matrix with collapsed expression data
#' 
#' @keywords internal
.collapse_replicates <- function(expression, groups) {
    unique_groups <- unique(groups)
    collapsed_list <- lapply(unique_groups, function(group) {
        group_samples <- groups == group
        if (sum(group_samples) == 1) {
            expression[, group_samples, drop = FALSE]
        } else {
            rowMeans(expression[, group_samples, drop = FALSE])
        }
    })
    
    result <- do.call(cbind, collapsed_list)
    colnames(result) <- as.character(unique_groups)
    return(result)
}

## METHOD IMPLEMENTATIONS

#' @title Convert BulkPhyloExpressionSet to Data Frame
#' @description Convert a BulkPhyloExpressionSet object back to the original data frame format
#' with phylostratum, gene ID, and expression data as columns.
#' 
#' @param phyex_set A BulkPhyloExpressionSet object
#' @param use_collapsed Logical indicating whether to use collapsed expression data (default: FALSE)
#' 
#' @return A data frame where column 1 contains phylostratum information, 
#' column 2 contains gene IDs, and columns 3+ contain expression data
#' 
#' @examples
#' # Convert BulkPhyloExpressionSet back to data frame
#' # df <- as_data_frame(bulk_phyex_set)
#' # df_collapsed <- as_data_frame(bulk_phyex_set, use_collapsed = TRUE)
#' 
#' @export
as_data_frame <- function(phyex_set, use_collapsed = FALSE) {
    if (use_collapsed) {
        expr_data <- phyex_set@expression_collapsed
    } else {
        expr_data <- phyex_set@expression
    }
    
    data.frame(
        Stratum = as.numeric(phyex_set@strata),
        GeneID = phyex_set@gene_ids,
        expr_data,
        check.names = FALSE
    )
}

#' @export
S7::method(collapse, BulkPhyloExpressionSet) <- function(phyex_set) {
    data <- tibble::tibble(Stratum=as.numeric(phyex_set@strata), 
                           GeneID=phyex_set@gene_ids, 
                           tibble::as_tibble(phyex_set@expression_collapsed))
    as_BulkPhyloExpressionSet(data)
}

#' @export
S7::method(select_genes, BulkPhyloExpressionSet) <- function(phyex_set, genes) {
    # Find indices of selected genes
    gene_indices <- match(genes, phyex_set@gene_ids)
    
    # Remove any NA indices (genes not found)
    valid_indices <- gene_indices[!is.na(gene_indices)]

    if (length(valid_indices) < length(gene_indices))
        warning("Some of the specified genes were not found in the dataset")
    
    if (length(valid_indices) == 0) {
        stop("None of the specified genes were found in the dataset")
    }
    
    # Create new data with selected genes
    data <- tibble::tibble(
        Stratum = phyex_set@strata[valid_indices],
        GeneID = phyex_set@gene_ids[valid_indices],
        tibble::as_tibble(phyex_set@expression[valid_indices, , drop = FALSE])
    )
    
    as_BulkPhyloExpressionSet(
        data = data,
        groups = phyex_set@groups,
        name = phyex_set@name,
        species = phyex_set@species,
        index_type = phyex_set@index_type
    )
}

#' @export
S7::method(print, BulkPhyloExpressionSet) <- function(x, ...) {
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
    
    # Print bulk-specific information
    cat("Number of samples:", x@num_samples, "\n")
    cat("Samples per condition:", table(x@groups), "\n")
}

# Aliases
as_PhyloExpressionSet <- as_BulkPhyloExpressionSet
