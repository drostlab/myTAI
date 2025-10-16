
#' @title Late Conservation Test
#' @description Test for late conservation patterns in transcriptomic data by comparing
#' late developmental stages to early and mid stages.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param modules A named list with elements 'early', 'mid', and 'late' containing
#' stage indices for each developmental module
#' @param ... Additional arguments passed to stat_generic_conservation_test
#' 
#' @return A ConservationTestResult object with late conservation test results
#' 
#' @details
#' The late conservation test evaluates whether later developmental stages show
#' lower transcriptomic index values (indicating older genes) compared to earlier
#' stages. The test computes a score based on the minimum difference between
#' early vs. late and mid vs. late TXI values.
#' 
#' @examples
#' # Define developmental modules
#' modules <- list(early = 1:2, mid = 3:5, late = 6:7)
#' result <- stat_late_conservation_test(example_phyex_set_old, modules)
#' 
#' @seealso \code{\link{stat_generic_conservation_test}}, \code{\link{stat_early_conservation_test}}
#' @export
stat_late_conservation_test <- function(phyex_set, modules, ...) {
    t <- stat_generic_conservation_test(phyex_set, 
                                        test_name="Late Conservation Test",
                                        scoring_function=\(x) lc_score(x, modules),
                                        fitting_dist=distributions$normal,
                                        alternative="greater",
                                        p_label="p_lct",
                                        ...)
    t@modules <- modules
    return(t)
}

#' @title Late Conservation Score Function
#' @description Compute the late conservation score by comparing early and mid
#' developmental stages to late stages.
#' 
#' @param txi Numeric vector of transcriptomic index values
#' @param modules A named list with elements 'early', 'mid', and 'late' containing
#' stage indices for each developmental module
#' 
#' @return A numeric value representing the late conservation score
#' 
#' @details
#' The score is computed as the minimum of:
#' - D1: mean(early) - mean(late)
#' - D2: mean(mid) - mean(late)
#' 
#' Higher scores indicate stronger late conservation patterns.
#' # Compute late conservation score
#' # modules <- list(early = 1:3, mid = 4:6, late = 7:9)
#' # score <- lc_score(txi_values, modules)
#' 
#' @keywords internal
lc_score <- function(txi, 
                     modules) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    D1 <- mean(txi[modules$early]) - mean(txi[modules$late])
    D2 <- mean(txi[modules$mid]) - mean(txi[modules$late])
    
    D_min <- min(D1, D2)
    
    return(D_min)
}