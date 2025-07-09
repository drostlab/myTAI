
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

ec_score <- function(txi, 
                     modules) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    D1 <- mean(txi[modules$mid]) - mean(txi[modules$early])
    D2 <- mean(txi[modules$late]) - mean(txi[modules$early])
    
    D_min <- min(D1, D2)
    
    return(D_min)
}



