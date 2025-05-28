
#' @import ggplot2
plot_gene_profiles <- function(phyex_set,
                               show_CI=T,
                               show_mean=TRUE,
                               colours=NULL) {
    
    if (is.null(colours))
        colours = rep("black", phyex_set@num_genes)
    
    df <- phyex_set@data |>
        tidyr::pivot_longer(-c(Stratum, GeneID), names_to="Sample", values_to="Expression") |>
        left_join(data.frame(Condition=factor(phyex_set@rep_groups, levels=unique(phyex_set@rep_groups)), Sample=phyex_set@sample_names), by="Sample") |>
        group_by(Stratum, GeneID, Condition) |>
        summarise(min=min(Expression), max=max(Expression), Expression=mean(Expression), .groups="drop")
    
    
    
    
    
    p <- ggplot(df, aes(x=factor(Condition, unique(Condition)),
                        group=GeneID,
                        colour=GeneID)) +
        geom_line(aes(y=Expression), alpha=0.5) +
        scale_color_manual(values=colours) +
        labs(x=phyex_set@conditions_label,
             y="Expression") +
        scale_x_discrete(labels = ~ stringr::str_wrap(., 10)) +
        theme_minimal()
    if (show_CI)
        p <- p + 
            geom_ribbon(aes(ymin=min, ymax=max), fill="gray", alpha=0.05)
    
    if (show_mean) {
        df_mean <- df |> group_by(Condition) |> summarise(Expression=mean(Expression), .groups="drop")
        p <- p + 
            geom_line(data = df_mean,
                      aes(y=Expression), colour="red", size=2, group=0)
    }
    
    
    return(p)
}