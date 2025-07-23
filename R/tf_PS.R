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
#' @examplesIf FALSE
#'
#' # get the relative expression profiles for each phylostratum
#' tfPES <- tf_PS(phyex_set, transform = "qr")
#' head(tfPES)
#'
#' @export

tfPS <- function(phyex_set, transform="qr") {
    if (!S7::S7_inherits(phyex_set, PhyloExpressionSetBase)) {
        stop("Input must be a PhyloExpressionSet S7 object.", call. = FALSE)
    }
    ps_vec <- as.numeric(phyex_set@strata)
    if (transform %in% c("qr", "quantilerank")) {
        ranks <- base::rank(ps_vec, ties.method = "average")
        tfPhylostratum <- (ranks - 0.5) / length(ps_vec)
    } else {
        stop("Choose the following transformation functions: \"qr\"")
    }
    
    # Get the original level names
    original_levels <- levels(phyex_set@strata)
    
    # Calculate the average transformed value for each unique phylostratum level
    unique_original_numeric <- sort(unique(ps_vec))
    level_transforms <- sapply(unique_original_numeric, function(level) {
        mean(tfPhylostratum[ps_vec == level])
    })
    
    # Create a mapping from original level names to transformed values
    level_mapping <- setNames(level_transforms, original_levels[unique_original_numeric])
    
    # Apply the transformation to each gene
    transformed_values <- level_mapping[as.character(phyex_set@strata)]
    
    # Create new factor with transformed values as levels and original names as labels
    unique_transformed <- sort(unique(transformed_values))
    original_names_sorted <- names(level_mapping)[order(level_mapping)]
    
    phyex_set@strata <- factor(transformed_values, 
                              levels = unique_transformed, 
                              labels = original_names_sorted)
    
    return(phyex_set)
}
