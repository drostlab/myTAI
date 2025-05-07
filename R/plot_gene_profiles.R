
#' @import ggplot2
plot_gene_profiles <- function(phyex_set,
                               show_CI=FALSE,
                               show_mean=TRUE,
                               group_by_strata=FALSE) {
    
    df <- phyex_set@data |>
        tidyr::pivot_longer(-c(Stratum, GeneID), names_to="Sample", values_to="Expression") |>
        left_join(data.frame(Condition=phyex_set@conditions, Sample=phyex_set@sample_names), by="Sample") |>
        group_by(Stratum, GeneID, Condition) |>
        summarise(min=min(Expression), Expression=mean(Expression), max=max(Expression), .groups="drop")
    
    p <- ggplot(df, aes(x=factor(Condition, unique(Condition)),
                        group=GeneID)) +
        geom_line(aes(y=Expression, colour=Stratum)) +
        labs(x=phyex_set@conditions_label,
             y="Expression") +
        scale_x_discrete(labels = ~ stringr::str_wrap(., 10)) +
        scale_colour_manual(values = PS_colours(phyex_set@num_strata)) +
        scale_fill_manual(values = PS_colours(phyex_set@num_strata)) +
        theme_minimal()
    if (show_CI)
        p <- p + 
            geom_ribbon(aes(ymin=min, ymax=max, fill=Stratum), alpha=0.05)
    
    if (show_mean) {
        df_mean <- df |> group_by(Condition) |> summarise(Expression=mean(Expression), .groups="drop")
        p <- p + 
            geom_line(data = df_mean,
                      aes(y=Expression), colour="red", size=2, group=0)
    }
        
    
    if (group_by_strata)
        p <- p + 
        facet_grid(. ~ Stratum) +
        labs(x=paste(phyex_set@conditions_label, "by Stratum")) +
        guides(x = "none", colour="none", fill="none")
    
    return(p)

}