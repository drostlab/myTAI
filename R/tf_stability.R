#' @title Perform Permutation Tests Under Different Transformations
#' @description \emph{tf_stability} statistically evaluates the
#' stability of phylotranscriptomics permutation tests (e.g., \code{ReductiveHourglassTest}, \code{FlatLineTest}, etc.)
#' under different data transformations using a \code{PhyloExpressionSet}.
#' @param phyex_set a \code{PhyloExpressionSet}.
#' @param conservation_test a conservation test function (e.g. \code{flatline_test}, \code{reductive_hourglass_test}, etc.)
#' @param transforms named list of transformation functions (default: \code{COUNT_TRANSFORMS})
#' @details
#' Assesses the stability of data transforms on the permutation test of choice.
#' See \code{\link{tf}}, \code{\link{ReductiveHourglassTest}}, \code{\link{FlatLineTest}}, etc.
#' @return Named numeric vector of p-values for each transformation.
#' @references
#' Lotharukpong JS et al. (2023) (unpublished)
#' @author Jaruwatana Sodai Lotharukpong
#' @export
tf_stability <- function(phyex_set,
                         conservation_test = flatline_test,
                         transforms = COUNT_TRANSFORMS) {
    # Validate input object
    if (!is.function(conservation_test)) {
        stop("conservation_test must be a function, e.g. flatline_test, reductive_hourglass_test, etc.", call. = FALSE)
    }

    vec_res <- numeric(length(transforms))
    names(vec_res) <- names(transforms)

    for (i in seq_along(transforms)) {
        tf_fun <- transforms[[i]]
        tf_name <- names(transforms)[i]
        tf_phyex <- transform_counts(phyex_set, FUN = tf_fun, FUN_name = tf_name)
        vec_res[i] <- conservation_test(tf_phyex, plot_result = FALSE)@p_value
    }
    return(vec_res)
}
