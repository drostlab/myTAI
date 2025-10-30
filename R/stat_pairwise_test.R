#' @title Pairwise Conservation Test
#' @description Test for significant differences in transcriptomic index values
#' between two groups of developmental stages.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param modules A named list with elements 'contrast1' and 'contrast2' containing
#' stage indices for each contrast group
#' @param alternative Character string specifying the alternative hypothesis:
#' "greater" (contrast1 > contrast2) or "less" (contrast1 < contrast2)
#' @param ... Additional arguments passed to stat_generic_conservation_test
#' 
#' @return A ConservationTestResult object with pairwise test results
#' 
#' @details
#' The pairwise test compares the mean transcriptomic index values between two
#' groups of developmental stages. This is useful for testing specific hypotheses
#' about differences in gene age composition between developmental periods.
#' 
#' @author Jaruwatana Sodai Lotharukpong
#' 
#' @examples
#' # Define contrast groups
#' modules <- list(contrast1 = 1:3, contrast2 = 4:7)
#' result <- stat_pairwise_test(example_phyex_set, modules, alternative = "greater")
#' 
#' @seealso \code{\link{stat_generic_conservation_test}}
#' @export
stat_pairwise_test <- function(phyex_set,
                          modules,
                          alternative = c("greater", "less"),
                          ...
) {
    alternative <- match.arg(alternative)
    t <- stat_generic_conservation_test(phyex_set, 
                                        test_name="Pairwise Test",
                                        scoring_function=\(x) pair_score(x, modules, alternative),
                                        fitting_dist=distributions$normal,
                                        alternative=alternative,
                                        p_label="p_pair",
                                        ...)
    t@modules <- modules
    return(t)
}

#' @title Pairwise Score Function
#' @description Compute the pairwise contrast score between two groups of developmental stages.
#' 
#' @param txi Numeric vector of transcriptomic index values
#' @param modules A named list with elements 'contrast1' and 'contrast2' containing
#' stage indices for each contrast group
#' @param alternative Character string specifying the alternative hypothesis:
#' "greater" or "less"
#' 
#' @return A numeric value representing the pairwise contrast score
#' 
#' @details
#' The score is computed as mean(contrast1) - mean(contrast2).
#' For alternative = "less", the score is negated.
#' # Compute pairwise score
#' # modules <- list(contrast1 = 1:3, contrast2 = 7:9)
#' # score <- pair_score(txi_values, modules, "greater")
#' 
#' @keywords internal
pair_score <- function(txi,
                       modules,
                       alternative = c("greater", "less")) {
    alternative <- match.arg(alternative)
    if(!setequal(c("contrast1", "contrast2"), names(modules)))
        stop('`modules` must have the structure: `list(contrast1 = ..., contrast2 = ...)`')
    D_constrast <- mean(txi[modules$contrast1]) - mean(txi[modules$contrast2])
    
    
    return(D_constrast)
}

