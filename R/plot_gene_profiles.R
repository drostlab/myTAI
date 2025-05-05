
#' @import ggplot2
plot_gene_profiles <- function(phyex_set,
                               show_CI=FALSE,
                               group_by_strata=FALSE) {
    
    df <- phyex_set@data |>
        tidyr::pivot_longer(-c(Stratum, GeneID), names_to="Sample", values_to="Expression") |>
        left_join(data.frame(Condition=phyex_set@conditions, Sample=phyex_set@sample_names), by="Sample") |>
        group_by(Stratum, GeneID, Condition) |>
        summarise(min=min(Expression), mean=mean(Expression), max=max(Expression), .groups="drop")
    
    p <- ggplot(df, aes(x=factor(Condition, unique(Condition)),
                        group=GeneID)) +
        geom_line(aes(y=mean, colour=Stratum)) +
        labs(x=phyex_set@conditions_label,
             y="Expression") +
        scale_x_discrete(labels = ~ stringr::str_wrap(., 10)) +
        scale_colour_manual(values = PS_colours(phyex_set@num_strata)) +
        scale_fill_manual(values = PS_colours(phyex_set@num_strata)) +
        theme_minimal()
    if (show_CI)
        p <- p + 
        geom_ribbon(aes(ymin=min, ymax=max, fill=Stratum), alpha=0.05)
        
    
    if (group_by_strata)
        p <- p + 
        facet_grid(. ~ Stratum) +
        guides(x = "none")
    
    return(p)

}