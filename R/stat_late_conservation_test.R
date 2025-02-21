
late_conservation_test <- function(phyex_set, modules) {
    return(generic_conservation_test(phyex_set, 
                                     test_name="Late Conservation Test",
                                     scoring_function=\(x) lc_score(x, modules),
                                     fitting_dist="normal",
                                     alternative="greater"))
}

lc_score <- function(txi, 
                     modules) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    D1 <- mean(txi[modules$early]) - mean(txi[modules$late])
    D2 <- mean(txi[modules$mid]) - mean(txi[modules$late])
    
    D_min <- min(D1, D2)
    
    return(D_min)
}