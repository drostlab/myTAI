
reverse_hourglass_test <- function(phyex_set, 
                                   modules) {
    return(generic_conservation_test(phyex_set, 
                                     test_name="Reverse Hourglass Test",
                                     scoring_function=\(x) reverse_hourglass_score(x, modules),
                                     fitting_dist="normal",
                                     alternative="greater"))
}

reverse_hourglass_score <- function(txi, 
                                    modules) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    D1 <- mean(txi[modules$mid]) - mean(txi[modules$early])
    D2 <- mean(txi[modules$mid]) - mean(txi[modules$late])
    
    D_min <- min(D1, D2)
    
    return(D_min)
}



