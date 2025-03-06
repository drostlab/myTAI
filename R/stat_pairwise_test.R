
pairwise_test <- function(phyex_set,
                          modules,
                          alternative = c("greater", "less"),
                          ...
) {
    alternative <- match.arg(alternative)
    t <- generic_conservation_test(phyex_set, 
                                   test_name="Pairwise Test",
                                   scoring_function=\(x) pair_score(x, modules, alternative),
                                   fitting_dist=distributions$normal,
                                   alternative=alternative,
                                   ...)
    t@modules <- modules
    return(t)
}

pair_score <- function(txi,
                       modules,
                       alternative = c("greater", "less")) {
    alternative <- match.arg(alternative)
    if(!setequal(c("contrast1", "contrast2"), names(modules)))
        stop('`modules` must have the structure: `list(contrast1 = ..., contrast2 = ...)`')
    D_constrast <- mean(txi[modules$contrast1]) - mean(txi[modules$contrast2])
    
    if(alternative == "less"){
        D_constrast <- -D_constrast
    }
    
    return(D_constrast)
}

