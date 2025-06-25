

plot_mean_var <- function(phyex_set) {
    df <- phyex_set@data_collapsed |>
        mutate(mean = rowMeans(phyex_set@counts_collapsed)) |>
        mutate(var = rowVars(phyex_set@counts_collapsed)) |>
        mutate(Stratum = as.numeric(Stratum))
    
    ggplot(df, aes(x=mean, y=var, colour=Stratum)) +
        geom_point(alpha=0.4) +
        scale_color_viridis_c() +
        theme_minimal()
}