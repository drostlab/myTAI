
#' @title Plot Phylostratum Contribution to Transcriptomic Index
#' @description Create a stacked area plot showing the contribution of each phylostratum
#' to the overall transcriptomic index across developmental stages or cell types.
#' 
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)
#' 
#' @return A ggplot2 object showing phylostratum contributions as a stacked area plot
#' 
#' @details
#' This function visualizes how different phylostrata contribute to the overall
#' transcriptomic index pattern across developmental stages (bulk data) or cell types
#' (single-cell data). Each area represents the contribution of a specific phylostratum, 
#' with older strata typically shown in darker colors and younger strata in lighter colors.
#' 
#' The plot uses the PS_colours function to create a consistent color scheme that
#' matches other myTAI visualizations.
#' 
#' @examples
#' # Create contribution plot for bulk data
#' # contrib_plot <- plot_contribution(bulk_phyex_set)
#' 
#' # Create contribution plot for single-cell data
#' # contrib_plot_sc <- plot_contribution(sc_phyex_set)
#' 
#' @import dplyr ggplot2
#' @export
plot_contribution <- function(phyex_set) {
    contribution <- sTXI(phyex_set, option="identity")
    strata_labels <- levels(phyex_set@strata)
    df <- data.frame(contribution) |>
        mutate(Stratum = factor(strata_labels, levels = rev(strata_labels))) |>
        tidyr::pivot_longer(-Stratum, names_to = "Identity", values_to = "TXI") |>
        mutate(Identity = factor(Identity, levels=unique(Identity)))
    
    
    p <- ggplot(df, 
                aes(x = Identity, 
                    y = TXI, 
                    fill = Stratum, 
                    group = Stratum)) +
        geom_area(position = "stack", colour="black", size=0.2) +
        labs(x=phyex_set@identities_label, y=phyex_set@index_full_name) +
        scale_fill_manual(values = rev(PS_colours(phyex_set@num_strata))) +
        guides(fill=guide_legend(override.aes = list(size = 1.0), keyheight = 0.5, keywidth = 0.5)) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(p)
}