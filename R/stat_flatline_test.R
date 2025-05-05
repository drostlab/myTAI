
flatline_test <- function(phyex_set, ...) {
    return(generic_conservation_test(phyex_set, 
                                     test_name="Flat Line Test",
                                     scoring_function=stats::var,
                                     fitting_dist=distributions$gamma,
                                     alternative="greater",
                                     p_label="p_flt",
                                     ...))
}