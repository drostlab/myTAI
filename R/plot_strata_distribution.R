
#' @import ggplot2 dplyr tidyr
plot_strata_distribution <- function(strata_vector,
                                     selected_gene_ids = names(strata_vector),
                                     as_log_obs_exp = FALSE
                                     ) {
    df <- data.frame(Stratum=strata_vector, GeneID=names(strata_vector))
    df_selected <- df |> filter(GeneID %in% selected_gene_ids)
    levels = sort(unique(strata_vector))
    if (!as_log_obs_exp) {
        ggplot(df_selected, aes(factor(Stratum,levels=levels), fill=factor(Stratum,levels=levels))) + 
            geom_bar() + 
            labs(fill="Stratum") +
            xlab("Stratum") + 
            ylab("Gene Count") +
            scale_x_discrete(drop = FALSE) +
            scale_fill_manual(values = PS_colours(length(levels))) +
            guides(fill="none") + 
            theme_minimal()
    }
    else {
        counts_all <- df |> count(Stratum) |> mutate(p_all = n / sum(n))
        counts_sel <- df_selected |> count(Stratum) |> mutate(p_sel = n / sum(n))
        log_ratio <- full_join(counts_all, counts_sel, by= "Stratum") |>
            mutate(log_obs_exp = log2((p_sel + 1e-6)/ (p_all + 1e-6)))
        
        ggplot(log_ratio, aes(factor(Stratum, levels=levels), 
                              y=log_obs_exp,
                              fill= log_obs_exp)
                              ) +
            geom_bar(stat = "identity") +
            labs(fill="Stratum") +
            xlab("Stratum") +
            ylab("Log(Obs/Exp)") +
            scale_x_discrete(drop = FALSE) +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
            guides(fill="none", colour="none") + 
            theme_minimal()
    
        
    }
    # TODO: integrate PlotEnrichment here
}