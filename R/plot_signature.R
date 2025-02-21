#' @import ggplot2
plot_signature <- function(phyex_set,
                           CI_low=.025,
                           CI_high=.975) {
    CIs <- TXI_conf_int(phyex_set, probs=c(CI_low, CI_high))
    df <- tibble::tibble(Stage = phyex_set@stages,
                         TXI = phyex_set@TXI,
                         CI_low=CIs[1, ],
                         CI_high=CIs[2, ]
    )
    p <- ggplot(data=df) +
        geom_line(
            aes(
                x=Stage,
                y=TXI,
                group=1),
            lwd=2
        ) +
        labs(
            x=phyex_set@stages_label,
            y=phyex_set@index_full_name
        )
    p <- p +
        geom_ribbon(
            aes(
                x=Stage,
                ymin=CI_low,
                ymax=CI_high,
                group=1
            ),
            fill="gray70",
            alpha=0.3
        )
    p <- p + theme_minimal()
    
    return(p)
}