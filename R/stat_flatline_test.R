#' @title Flat Line Test for Conservation Pattern
#' @description Perform a flat line test to assess whether the transcriptomic index 
#' profile shows a flat (non-varying) pattern across developmental stages.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param ... Additional arguments passed to stat_generic_conservation_test
#' 
#' @return A test result object containing p-value and test statistics
#' 
#' @details
#' The flat line test evaluates whether the TXI profile remains constant across
#' developmental stages by testing the variance of the profile against a null
#' distribution. A significant result indicates rejection of the flat line pattern.
#' 
#' @examples
#' # Perform flat line test
#' result <- stat_flatline_test(example_phyex_set)
#' 
#' @seealso \code{\link{stat_generic_conservation_test}}
#' @export
stat_flatline_test <- function(phyex_set, ...) {
    check_PhyloExpressionSet(phyex_set)
    return(stat_generic_conservation_test(phyex_set, 
                                     test_name="Flat Line Test",
                                     scoring_function=stats::var,
                                     fitting_dist=distributions$gamma,
                                     alternative="greater",
                                     p_label="p_flt",
                                     ...))
}