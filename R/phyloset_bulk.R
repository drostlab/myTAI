#' @title Bulk PhyloExpressionSet Class
#' @description S7 class for bulk RNA-seq phylotranscriptomic expression data.
#' This class handles expression data with biological replicates and provides
#' bootstrapping functionality for statistical analysis.
#' 
#' @param strata Factor vector of phylostratum assignments for each gene
#' @param strata_values Numeric vector of phylostratum values used in TXI calculations
#' @param expression Matrix of expression counts with genes as rows and samples as columns
#' @param groups Factor vector indicating which identity each sample belongs to
#' @param name Character string naming the dataset (default: "Phylo Expression Set")
#' @param species Character string specifying the species (default: NULL)
#' @param index_type Character string specifying the transcriptomic index type (default: "TXI")
#' @param identities_label Character string labeling the identities (default: "Stages")
#' @param null_conservation_sample_size Numeric value for null conservation sample size (default: 5000)
#' @param .null_conservation_txis Precomputed null conservation TXI values (default: NULL)
#' @param .bootstrapped_txis Precomputed bootstrapped TXI values (default: NULL)
#' 
#' @return A BulkPhyloExpressionSet object
#' 
#' @details
#' The BulkPhyloExpressionSet class is designed for bulk RNA-seq data with biological replicates.
#' It extends the base PhyloExpressionSetBase class with bulk-specific functionality.
#' 
#' \strong{Replicate Handling:}
#' Expression data across biological replicates is collapsed by taking row means within each
#' experimental condition or developmental stage.
#' 
#' \strong{Computed Properties:}
#' In addition to inherited computed properties from the base class, this class provides:
#' \itemize{
#'   \item \code{expression_collapsed} - Matrix of expression data collapsed across replicates (genes x identities)
#'   \item \code{bootstrapped_txis} - Matrix of bootstrapped TXI values for statistical inference (500 bootstrap samples x identities)
#' }
#' 
#' Inherited computed properties from PhyloExpressionSetBase include:
#' \itemize{
#'   \item \code{gene_ids} - Character vector of gene identifiers
#'   \item \code{identities} - Character vector of identity labels  
#'   \item \code{sample_names} - Character vector of sample names
#'   \item \code{num_identities} - Integer count of unique identities
#'   \item \code{num_samples} - Integer count of total samples
#'   \item \code{num_genes} - Integer count of genes
#'   \item \code{num_strata} - Integer count of phylostrata
#'   \item \code{index_full_name} - Full name of the transcriptomic index type
#'   \item \code{group_map} - List mapping identity names to sample names
#'   \item \code{TXI} - Numeric vector of TXI values for each identity
#'   \item \code{TXI_sample} - Numeric vector of TXI values for each sample
#'   \item \code{null_conservation_txis} - Matrix of null conservation TXI values for statistical testing
#' }
#' 
#' \strong{Statistical Analysis:}
#' The class supports confidence interval estimation and standard deviation calculation
#' through bootstrapped TXI values, enabling robust statistical analysis of developmental
#' or experimental patterns.
#' 
#' @import S7
#' @export
BulkPhyloExpressionSet <- new_class("BulkPhyloExpressionSet",
    parent = PhyloExpressionSetBase,
    properties = list(
        # Implemented abstract properties
        expression_collapsed = new_property(
            getter = function(self) .collapse_replicates(self@expression, self@groups)
        ),
        identities_label = new_property(
            class = class_character,
            default = "Stages"
        ),
        # Bulk specific properties
        .bootstrapped_txis = new_property(
            default = NULL
        ),
        bootstrapped_txis = new_property(
            getter = function(self) {
                if (is.null(self@.bootstrapped_txis)) {
                    # Compute and cache the result using expression_collapsed
                    ptxi <- .pTXI(self@expression_collapsed, self@strata_values)
                    computed_boot <- memo_generate_bootstrapped_txis(
                        ptxi,
                        self@expression_collapsed,
                        500
                    )
                    self@.bootstrapped_txis <- computed_boot
                    return(computed_boot)
                } else {
                    return(self@.bootstrapped_txis)
                }
            }
        )
    ),
    validator = function(self) {
        # This validation should/could be done in PhyloExpressionSetBase, but doesn't work because of a bug: 
        # https://github.com/RConsortium/S7/issues/539

        # Validate expression_collapsed rownames match gene_ids
        if (!identical(rownames(self@expression_collapsed), self@gene_ids)) 
            return("@expression_collapsed rownames must match @gene_ids")
        
        unique_groups <- sort(as.character(unique(self@groups)))
        if (!identical(sort(self@identities), unique_groups)) {
            return("@identities must match unique values in @groups")
        }
    }
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
#' 
#' @export
BulkPhyloExpressionSet_from_df <- function(data, 
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
    strata_values <- as.numeric(data[[1]])

    names(strata) <- gene_ids
    names(strata_values) <- gene_ids

    groups <- factor(groups, levels=unique(groups))
    
    expression <- as.matrix(data[3:ncol(data)])
    rownames(expression) <- gene_ids
    
    obj <- BulkPhyloExpressionSet(
        strata = strata,
        strata_values = strata_values,
        expression = expression,
        groups = groups,
        name = name,
        ...
    )


    return(obj)
}

#' @title as_BulkPhyloExpressionSet
#' @description
#' This function is an alias for \code{BulkPhyloExpressionSet_from_df}. 
#' Please refer to the documentation for \code{BulkPhyloExpressionSet_from_df} for usage details, arguments, and examples.
#' @param data A data frame with phylostratum assignments and gene expression data
#' @param groups Vector of group labels for the samples/replicates
#' @param name Character string to name the dataset
#' @param strata_legend Optional data frame mapping phylostratum numbers to labels
#' @param ... Additional arguments passed to BulkPhyloExpressionSet_from_df
#' @seealso \code{\link{BulkPhyloExpressionSet_from_df}}
#' @return A BulkPhyloExpressionSet object
#' @export
as_BulkPhyloExpressionSet <- BulkPhyloExpressionSet_from_df

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
                      groups = colnames(data[, 2:ncol(data)]),
                      name = NULL,
                      ...) {
    if (is.null(name)) name <- deparse(substitute(data))
    colnames(phylomap) <- c("Stratum", "GeneID")
    
    full_data <- data |>
        inner_join(phylomap, by="GeneID") |>
        relocate(Stratum, .before = 1)
    return(as_BulkPhyloExpressionSet(full_data, groups = groups, name = name, ...))
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
#' df <- as_data_frame(example_phyex_set)
#' df_collapsed <- as_data_frame(example_phyex_set, use_collapsed = TRUE)
#' 
#' @export
as_data_frame <- function(phyex_set, use_collapsed = FALSE) {
    if (use_collapsed) {
        expr_data <- phyex_set@expression_collapsed
    } else {
        expr_data <- phyex_set@expression
    }
    
    data.frame(
        Stratum = phyex_set@strata_values,
        GeneID = phyex_set@gene_ids,
        expr_data,
        check.names = FALSE
    )
}

#' @export
S7::method(collapse, BulkPhyloExpressionSet) <- function(phyex_set) {
    obj <- S7::valid_eventually(phyex_set, function(x) {
        ec <- x@expression_collapsed
        x@groups <- factor(colnames(ec), 
                           levels=colnames(ec))
        x@expression <- ec
        x
    })
    obj
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

    obj <- S7::valid_eventually(phyex_set, function(x) {
        x@strata <- phyex_set@strata[valid_indices]
        x@strata_values <- phyex_set@strata_values[valid_indices]
        x@expression <- phyex_set@expression[valid_indices, , drop = FALSE]
        x
    })
    
    obj
}

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

#' @title Check if object is a BulkPhyloExpressionSet
#' @description Checks if the input is a PhyloExpressionSet S7 object and throws an error if not.
#' @param phyex_set An object to check
#' @return Invisibly returns TRUE if check passes, otherwise throws an error
#' @export
check_BulkPhyloExpressionSet <- function(phyex_set) {
    if (!S7::S7_inherits(phyex_set, BulkPhyloExpressionSet)) {
        stop("Input must be a BulkPhyloExpressionSet S7 object.", call. = FALSE)
    }
    invisible(TRUE)
}

#' @title Confidence Intervals for Transcriptomic Index (TXI)
#' @description Compute confidence intervals for the TXI using bootstrapped TXI values.
#' @param phyex_set A BulkPhyloExpressionSet object
#' @param probs Numeric vector of probabilities for the confidence interval (default: c(0.025, 0.975))
#' @return A tibble with first column Identity names, second column lower bound, third column upper bound
#' @details
#' This function returns confidence intervals for the TXI for each identity (sample or group),
#' based on the bootstrapped TXI values stored in the PhyloExpressionSet object.
#' @export
TXI_conf_int <- function(phyex_set, 
                         probs=c(.025, .975)) {
    check_BulkPhyloExpressionSet(phyex_set)
    CIs <- apply(phyex_set@bootstrapped_txis, 2, quantile, probs=probs)
    tibble::tibble(
        Identity = factor(phyex_set@identities, levels = unique(as.character(phyex_set@identities))),
        lb = as.numeric(CIs[1,]),
        ub = as.numeric(CIs[2,])
    )
}

#' @title Standard Deviation for TXI
#' @description Return a named vector of standard deviations for the TXI for each identity.
#' @param phyex_set A BulkPhyloExpressionSet object
#' @return Named numeric vector of standard deviations, names are identities
#' @export
TXI_std_dev <- function(phyex_set) {
    check_BulkPhyloExpressionSet(phyex_set)
    std_dev <- apply(phyex_set@bootstrapped_txis, 2, sd)
    names(std_dev) <- as.character(phyex_set@identities)
    std_dev
}