#' @title Plot Mean-Variance Relationship
#' @description Create a scatter plot showing the relationship between mean expression
#' and variance for genes, colored by phylostratum.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return A ggplot2 object showing the mean-variance relationship
#' 
#' @details
#' This function plots the mean expression versus variance for each gene,
#' with points colored by phylostratum. This visualization helps identify
#' expression patterns and heteroscedasticity in the data.
#' 
#' @examples
#' # Create mean-variance plot
#' # mv_plot <- plot_mean_var(phyex_set)
#' 
#' @import ggplot2
#' @import dplyr
#' @export
plot_mean_var <- function(phyex_set) {
    df <- phyex_set@data_collapsed |>
        mutate(mean = rowMeans(phyex_set@counts_collapsed)) |>
        mutate(var = rowVars(phyex_set@counts_collapsed))
    
    ggplot(df, aes(x=mean, y=var, colour=Stratum)) +
        geom_point(alpha=0.7, size=0.2) +
        scale_color_manual(values = PS_colours(length(unique(df$Stratum))), na.value = "grey50", name = "Stratum") +
        theme_minimal()
}