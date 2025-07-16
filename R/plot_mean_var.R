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