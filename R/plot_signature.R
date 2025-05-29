
#' @import ggplot2
plot_signature <- function(phyex_set,
                           show_CI=FALSE,
                           show_bootstraps=FALSE,
                           CI_low=.025,
                           CI_high=.975,
                           show_reps=TRUE,
                           show_p_val=FALSE,
                           conservation_test=flatline_test,
                           colour = NULL,
                           ...) {
    p <- ggplot()
    
    # Plot CI (unless showing bootstraps)
    if (show_CI && !show_bootstraps) {
        CI <- TXI_conf_int(phyex_set, low_q=CI_low, high_q=CI_high)
        geom <- if (phyex_set@is_time_series) geom_ribbon else geom_linerange
        p <- p + geom(data=tibble::tibble(Condition = phyex_set@conditions,
                                                 CI_low=CI$low, 
                                                 CI_high=CI$high),
                             aes(x=Condition, 
                                 ymin=CI_low, 
                                 ymax=CI_high, 
                                 group=0,
                                 fill=phyex_set@name)
                             , alpha=0.3)
    }
    
    # Plot bootstraps
    if (show_bootstraps) {
        df_sample <- as_tibble(phyex_set@bootstrapped_txis) |>
            rowid_to_column("Id") |>
            tidyr::pivot_longer(cols=phyex_set@conditions, names_to="Condition", values_to = "TXI")
        geom <- if (phyex_set@is_time_series) geom_line else purrr::partial(geom_jitter, width=0.05)
        p <- p + geom(data=df_sample,
                           aes(x=Condition, 
                               y=TXI, 
                               group=Id,
                               colour=phyex_set@name),
                           alpha=0.04)
    }
    
    # Plot TXI
    df <- tibble::tibble(Condition = phyex_set@conditions,
                         TXI = phyex_set@TXI)
    if (phyex_set@is_time_series) {
        p <- p + 
            geom_line(data=df, 
                      aes(x=Condition, 
                          y=TXI, 
                          group=0), 
                      colour="black",
                      lwd=2.3,
                      lineend="round") +
            geom_line(data=df, 
                      aes(x=Condition, 
                          y=TXI, 
                          group=0,
                          colour=phyex_set@name), 
                      lwd=1.5,
                      lineend="round") 
    }
    else {
        p <- p + geom_point(data=df,
                            aes(x=Condition, 
                                y=TXI,
                                fill=phyex_set@name),
                            shape=21, colour="black", size=2, stroke=0.8)
    }
    
    p <- p + 
        labs(
            x=phyex_set@conditions_label,
            y=phyex_set@index_full_name
        ) +
        guides(colour="none", fill="none") +
        theme_minimal()
    
    # Plot replicates
    if (show_reps) {
        df <- tibble::tibble(Condition = factor(phyex_set@rep_groups, levels=unique(phyex_set@rep_groups)),
                             TXI = phyex_set@TXI_reps)
        p <- p + geom_jitter(data=df,
                             aes(x=Condition, 
                                 y=TXI,
                                 fill=phyex_set@name),
                             shape=21, colour="black", size=2, stroke=0.8, width=0.05)
    }
    
    # Show p value
    if (show_p_val) {
        t <- conservation_test(phyex_set, plot_result=FALSE)
        label <- paste(t@p_label, "=", signif(t@p_value, 3))
        p <- p + 
            annotate("label", label=label, fill="white",
                     x=phyex_set@num_conditions * 0.7, y=mean(phyex_set@TXI_reps) + 0.1)
    }
    
    if (!is.null(colour)) {
        p <- p +
            scale_color_manual(values = c(colour)) +
            scale_fill_manual(values = c(colour))
    }
    
    return(p)
}





