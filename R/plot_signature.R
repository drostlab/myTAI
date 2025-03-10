

plot_signature <- S7::new_generic("plot_signature", "phyex_set")

#TODO display p value as label
#' @import ggplot2
S7::method(plot_signature, PhyloExpressionSet) <- function(phyex_set,
                                                           show_CI=TRUE,
                                                           show_bootstraps=FALSE,
                                                           CI_probs=c(.025, .975)) {
    p <- ggplot()
    
    # Plot CI
    if (show_CI) {
        CIs <- TXI_conf_int(phyex_set, probs=CI_probs)
        p <- p + geom_ribbon(data=tibble::tibble(Condition = phyex_set@conditions,
                                                 CI_low=CIs[1, ], 
                                                 CI_high=CIs[2, ]),
                             aes(x=Condition, ymin=CI_low, ymax=CI_high, group=0),
                             fill="gray70", alpha=0.3)
    }
    
    # Plot bootstraps
    if (show_bootstraps) {
        df_sample <- as_tibble(phyex_set@bootstrapped_txis) |>
            rowid_to_column("Id") |>
            tidyr::pivot_longer(cols=phyex_set@conditions, names_to="Condition", values_to = "TXI")
        p <- p + geom_line(data=df_sample,
                           aes(x=Condition, y=TXI, group=Id),
                           alpha=0.2, colour="gray70")
    }
    
    # Plot TXI
    df <- tibble::tibble(Condition = phyex_set@conditions,
                         TXI = phyex_set@TXI)
    p <- p + geom_line(data=df, aes(x=Condition, y=TXI, group=0), lwd=2)
    
    p <- p + 
        labs(
            x=phyex_set@conditions_label,
            y=phyex_set@index_full_name
        ) +
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
        p <- p + geom_point(data=df,
                            aes(x=Condition, y=TXI))
    }
    return(p)
}

#' @import ggplot2
S7::method(plot_signature, PhyloExpressionSetMasked) <- function(phyex_set,
                                                                     ...) {
    
    # Plot the original pattern
    p <- plot_signature(phyex_set@full_set, ...)
    
    # Plot the perturbed pattern
    df <- tibble::tibble(Condition = phyex_set@conditions,
                         TXI = phyex_set@TXI)
    p <- p + geom_line(data=df, aes(x=Condition, y=TXI, group=0), lwd=1, colour="red")
    
    return(p)
}



