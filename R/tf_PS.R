#' @title Transform Phylostratum Values
#' @description
#' This function performs transformation of phylostratum values.
#'
#' @param phyex_set a \code{PhyloExpressionSet} S7 object.
#' @param transform a character vector of any valid function that transforms PS values.
#' Possible values can be:
#' \itemize{
#' \item \code{transform} = \code{"qr"} (or \code{"quantilerank"}) :
#' quantile rank transformation analogous to Julia function \code{StatsBase.quantilerank}
#' using \code{method = :tied}.
#' }
#' @details This function transforms the phylostratum assignment.
#' The return value of this function is a PhyloExpressionSet object with
#' transformed phylostratum \code{tfPhylostratum} as the first column, satisfying
#' \code{\link{PhyloExpressionSetBase}}. Note that the input \code{transform} must be an
#' available function, currently limited to only \code{"qr"} (or \code{"quantilerank"}).
#' @return a \code{PhyloExpressionSet} object storing transformed Phylostratum levels.
#' @author Jaruwatana Sodai Lotharukpong and Lukas Maischak
#' @seealso \code{\link{tf}}
#' @examples
#'
#' # get the relative expression profiles for each phylostratum
#' tfPES <- tf_PS(example_phyex_set, transform = "qr")
#'
#' @export

tf_PS <- function(phyex_set, transform="qr") {
    if (!S7::S7_inherits(phyex_set, PhyloExpressionSetBase)) {
        stop("Input must be a PhyloExpressionSet S7 object.", call. = FALSE)
    }
    ps_vec <- phyex_set@strata_values
    if (transform %in% c("qr", "quantilerank")) {
        ranks <- base::rank(ps_vec, ties.method = "average")
        tfPhylostratum <- (ranks - 0.5) / length(ps_vec)
    } else {
        stop("Choose the following transformation functions: \"qr\"")
    }
    
    # Get the original strata information
    original_strata <- phyex_set@strata
    original_strata_values <- phyex_set@strata_values
    original_levels <- levels(original_strata)
    
    # Calculate the average transformed value for each unique phylostratum level
    unique_original_values <- sort(unique(original_strata_values))
    level_transforms <- sapply(unique_original_values, function(level) {
        mean(tfPhylostratum[original_strata_values == level])
    })
    
    # Create a new factor with transformed values as levels but original labels preserved
    # The key insight: we need to map each gene's transformed value to the appropriate label
    
    # For each gene, get its transformed value
    gene_transformed_values <- sapply(original_strata_values, function(val) {
        level_transforms[which(unique_original_values == val)]
    })
    
    # Create mapping from unique original values to their labels
    # Since levels are in order 1, 2, 3, ..., we can map directly
    value_to_label <- setNames(original_levels, seq_len(length(original_levels)))
    
    # Get the label for each gene based on its original strata value
    gene_labels <- value_to_label[as.character(original_strata_values)]
    
    # Create sorted unique transformed values and corresponding labels
    sorted_indices <- order(unique_original_values)
    unique_transformed_sorted <- level_transforms[sorted_indices]
    labels_sorted <- original_levels[unique_original_values[sorted_indices]]
    
    # Create the new factor
    phyex_set@strata <- factor(gene_transformed_values, 
                              levels = unique_transformed_sorted, 
                              labels = labels_sorted)
    
    # Set the transformed values directly for each gene
    phyex_set@strata_values <- tfPhylostratum
    names(phyex_set@strata_values) <- phyex_set@gene_ids
    
    return(phyex_set)
}
