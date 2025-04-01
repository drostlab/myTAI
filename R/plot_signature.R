

plot_signature <- S7::new_generic("plot_signature", "phyex_set")

#TODO display p value as label
#' @import ggplot2
S7::method(plot_signature, PhyloExpressionSet) <- function(phyex_set,
                                                           show_CI=TRUE,
                                                           show_bootstraps=FALSE,
                                                           CI_probs=c(.025, .975),
                                                           ...) {
    p <- ggplot()
    
    # Plot CI (unless showing bootstraps)
    if (show_CI && !show_bootstraps) {
        CIs <- TXI_conf_int(phyex_set, probs=CI_probs)
        geom <- if (phyex_set@is_time_series) geom_ribbon else geom_linerange
        p <- p + geom(data=tibble::tibble(Condition = phyex_set@conditions,
                                                 CI_low=CIs[1, ], 
                                                 CI_high=CIs[2, ]),
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
    
    return(p)
}

#' @import ggplot2
S7::method(plot_signature, PhyloExpressionSetReplicates) <- function(phyex_set,
                                                                     ...,
                                                                     show_replicates=TRUE) {
    p <- plot_signature(S7::super(phyex_set, to=PhyloExpressionSet), ...)
    
    # Plot replicates
    if (show_replicates) {
        df <- tibble::tibble(Condition = factor(phyex_set@groups, levels=unique(phyex_set@groups)),
                             TXI = phyex_set@TXI_reps)
        p <- p + geom_jitter(data=df,
                            aes(x=Condition, 
                                y=TXI,
                                fill=phyex_set@name),
                            shape=21, colour="black", size=2, stroke=0.8, width=0.05)
    }
    return(p)
}

#' @import ggplot2
S7::method(plot_signature, PhyloExpressionSetMasked) <- function(phyex_set,
                                                                 show_full_set = TRUE,
                                                                     ...) {
    full_set <- phyex_set@full_set
    masked_set <- S7::convert(phyex_set, PhyloExpressionSet)
    # Plot the original pattern
    if (show_full_set)
        p <- plot_signature_multiple(c(full_set, masked_set))
    
    else
        p <- plot_signature(masked_set)
        
    return(p)
}



