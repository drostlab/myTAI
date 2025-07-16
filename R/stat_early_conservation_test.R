
#' @title Early Conservation Test
#' @description Test for early conservation patterns in transcriptomic data by comparing
#' early developmental stages to mid and late stages.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param modules A named list with elements 'early', 'mid', and 'late' containing
#' stage indices for each developmental module
#' @param ... Additional arguments passed to generic_conservation_test
#' 
#' @return A ConservationTestResult object with early conservation test results
#' 
#' @details
#' The early conservation test evaluates whether early developmental stages show
#' lower transcriptomic index values (indicating older genes) compared to later
#' stages. The test computes a score based on the minimum difference between
#' mid vs. early and late vs. early TXI values.
#' 
#' @examples
#' # Define developmental modules
#' # modules <- list(early = 1:3, mid = 4:6, late = 7:9)
#' # result <- early_conservation_test(phyex_set, modules)
#' 
#' @seealso \code{\link{generic_conservation_test}}, \code{\link{late_conservation_test}}
#' @export
early_conservation_test <- function(phyex_set, modules, ...) {
    t <- generic_conservation_test(phyex_set, 
                                   test_name="Early Conservation Test",
                                   scoring_function=\(x) ec_score(x, modules),
                                   fitting_dist=distributions$normal,
                                   alternative="greater",
                                   p_label = "p_ect",
                                   ...)
    t@modules <- modules
    return(t)
}

#' @title Early Conservation Score Function
#' @description Compute the early conservation score by comparing mid and late
#' developmental stages to early stages.
#' 
#' @param txi Numeric vector of transcriptomic index values
#' @param modules A named list with elements 'early', 'mid', and 'late' containing
#' stage indices for each developmental module
#' 
#' @return A numeric value representing the early conservation score
#' 
#' @details
#' The score is computed as the minimum of:
#' - D1: mean(mid) - mean(early)
#' - D2: mean(late) - mean(early)
#' 
#' Higher scores indicate stronger early conservation patterns.
#' 
#' @examples
#' # Compute early conservation score
#' # modules <- list(early = 1:3, mid = 4:6, late = 7:9)
#' # score <- ec_score(txi_values, modules)
#' 
#' @keywords internal
ec_score <- function(txi, 
                     modules) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    D1 <- mean(txi[modules$mid]) - mean(txi[modules$early])
    D2 <- mean(txi[modules$late]) - mean(txi[modules$early])
    
    D_min <- min(D1, D2)
    
    return(D_min)
}



