
reductive_hourglass_test <- function(phyex_set, modules, ...) {
    t <- generic_conservation_test(phyex_set, 
                                   test_name="Reductive Hourglass Test",
                                   scoring_function=\(x) reductive_hourglass_score(x, modules),
                                   fitting_dist=distributions$normal,
                                   alternative="greater",
                                   ...)
    t@modules <- modules
    return(t)
}

reductive_hourglass_score <- function(txi, 
                                      modules) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    D1 <- mean(txi[modules$early]) - mean(txi[modules$mid])
    D2 <- mean(txi[modules$late]) - mean(txi[modules$mid])
    
    D_min <- min(D1, D2)
    
    return(D_min)
}



