
#' @import dplyr ggplot2
#' @export
plot_contribution <- function(phyex_set) {
    contribution <- sTXI(phyex_set, option="identity")
    df <- data.frame(contribution) |>
        mutate(PS = factor(row_number())) |>
        tidyr::pivot_longer(-PS, names_to = "Condition", values_to = "TXI") |>
        mutate(Condition = factor(Condition, levels=unique(Condition))) |>
        mutate(PS = forcats::fct_rev(PS)) # reverse PS order
    
    
    p <- ggplot(df, 
                aes(x = Condition, 
                    y = TXI, 
                    fill = PS, 
                    group = PS)) +
        geom_area(position = "stack", colour="black", size=0.2) +
        labs(x=phyex_set@conditions_label, y=phyex_set@index_full_name) +
        scale_fill_manual(values = rev(PS_colours(phyex_set@num_strata))) +
        guides(fill=guide_legend(override.aes = list(size = 1.0), keyheight = 0.5, keywidth = 0.5)) +
        theme_minimal()
    return(p)
}