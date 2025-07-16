
#' @title Reductive Hourglass Test
#' @description Test for reductive hourglass patterns in transcriptomic data by comparing
#' early and late developmental stages to mid developmental stages.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param modules A named list with elements 'early', 'mid', and 'late' containing
#' stage indices for each developmental module
#' @param ... Additional arguments passed to generic_conservation_test
#' 
#' @return A ConservationTestResult object with reductive hourglass test results
#' 
#' @details
#' The reductive hourglass test evaluates whether mid developmental stages show
#' lower transcriptomic index values (indicating older genes) compared to both
#' early and late stages. This creates an hourglass-shaped pattern where ancient
#' genes dominate during mid-development. The test computes a score based on the
#' minimum difference between early vs. mid and late vs. mid TXI values.
#' 
#' @examples
#' # Define developmental modules
#' # modules <- list(early = 1:3, mid = 4:6, late = 7:9)
#' # result <- reductive_hourglass_test(phyex_set, modules)
#' 
#' @seealso \code{\link{generic_conservation_test}}, \code{\link{reverse_hourglass_test}}
#' @export
reductive_hourglass_test <- function(phyex_set, modules, ...) {
    t <- generic_conservation_test(phyex_set, 
                                   test_name="Reductive Hourglass Test",
                                   scoring_function=\(x) reductive_hourglass_score(x, modules),
                                   fitting_dist=distributions$normal,
                                   alternative="greater",
                                   p_label="p_red",
                                   ...)
    t@modules <- modules
    return(t)
}

#' @title Reductive Hourglass Score Function
#' @description Compute the reductive hourglass score by comparing early and late
#' developmental stages to mid stages.
#' 
#' @param txi Numeric vector of transcriptomic index values
#' @param modules A named list with elements 'early', 'mid', and 'late' containing
#' stage indices for each developmental module
#' 
#' @return A numeric value representing the reductive hourglass score
#' 
#' @details
#' The score is computed as the minimum of:
#' - D1: mean(early) - mean(mid)
#' - D2: mean(late) - mean(mid)
#' 
#' Higher scores indicate stronger reductive hourglass patterns (mid stages
#' dominated by older genes).
#' 
#' @examples
#' # Compute reductive hourglass score
#' # modules <- list(early = 1:3, mid = 4:6, late = 7:9)
#' # score <- reductive_hourglass_score(txi_values, modules)
#' 
#' @keywords internal
reductive_hourglass_score <- function(txi, 
                                      modules) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    D1 <- mean(txi[modules$early]) - mean(txi[modules$mid])
    D2 <- mean(txi[modules$late]) - mean(txi[modules$mid])
    
    D_min <- min(D1, D2)
    
    return(D_min)
}



