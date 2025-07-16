
#' @title Reverse Hourglass Test
#' @description Test for reverse hourglass patterns in transcriptomic data by comparing
#' mid developmental stages to early and late developmental stages.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param modules A named list with elements 'early', 'mid', and 'late' containing
#' stage indices for each developmental module
#' @param ... Additional arguments passed to generic_conservation_test
#' 
#' @return A ConservationTestResult object with reverse hourglass test results
#' 
#' @details
#' The reverse hourglass test evaluates whether mid developmental stages show
#' higher transcriptomic index values (indicating younger genes) compared to both
#' early and late stages. This creates a reverse hourglass pattern where recently
#' evolved genes dominate during mid-development. The test computes a score based
#' on the minimum difference between mid vs. early and mid vs. late TXI values.
#' 
#' @examples
#' # Define developmental modules
#' # modules <- list(early = 1:3, mid = 4:6, late = 7:9)
#' # result <- reverse_hourglass_test(phyex_set, modules)
#' 
#' @seealso \code{\link{generic_conservation_test}}, \code{\link{reductive_hourglass_test}}
#' @export
reverse_hourglass_test <- function(phyex_set, 
                                   modules,
                                   ...) {
    t <- generic_conservation_test(phyex_set, 
                                   test_name="Reverse Hourglass Test",
                                   scoring_function=\(x) reverse_hourglass_score(x, modules),
                                   fitting_dist=distributions$normal,
                                   alternative="greater",
                                   p_label="p_rht",
                                   ...)
    t@modules <- modules
    return(t)
}

#' @title Reverse Hourglass Score Function
#' @description Compute the reverse hourglass score by comparing mid developmental
#' stages to early and late stages.
#' 
#' @param txi Numeric vector of transcriptomic index values
#' @param modules A named list with elements 'early', 'mid', and 'late' containing
#' stage indices for each developmental module
#' 
#' @return A numeric value representing the reverse hourglass score
#' 
#' @details
#' The score is computed as the minimum of:
#' - D1: mean(mid) - mean(early)
#' - D2: mean(mid) - mean(late)
#' 
#' Higher scores indicate stronger reverse hourglass patterns (mid stages
#' dominated by younger genes).
#' 
#' @examples
#' # Compute reverse hourglass score
#' # modules <- list(early = 1:3, mid = 4:6, late = 7:9)
#' # score <- reverse_hourglass_score(txi_values, modules)
#' 
#' @keywords internal
reverse_hourglass_score <- function(txi, 
                                    modules) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    D1 <- mean(txi[modules$mid]) - mean(txi[modules$early])
    D2 <- mean(txi[modules$mid]) - mean(txi[modules$late])
    
    D_min <- min(D1, D2)
    
    return(D_min)
}



