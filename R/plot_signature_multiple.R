


#' @import ggplot2 purrr
plot_signature_multiple <- function(phyex_sets,
                                    legend_title="Phylo Expression Set",
                                    ...) {
    
    layers <- phyex_sets |>
        map(plot_signature, ...) |>
        map(~ .x$layers) |>
        unlist(c()) |>
        .sort_layers()
    
    p <- ggplot() +
        layers +
        labs(
            x=phyex_sets[[1]]@conditions_label,
            y=phyex_sets[[1]]@index_full_name,
            colour=legend_title,
            fill=legend_title
        ) +
        theme_minimal()
        
    return(p)
}

# a hacky way to sort the layers so that the TXI line is shown above
# the CI and bootstraps
.sort_layers <- function(layers) {
    is_round_line <- function(l) {
        inherits(l$geom, "GeomLine") && "lineend" %in% names(l$geom_params) && l$geom_params$lineend == "round"
    }
    
    is_point <- function(l) inherits(l$geom, "GeomPoint")
    
    others <- layers[!sapply(layers, function(l) is_round_line(l) || is_point(l))]
    round_lines <- layers[sapply(layers, is_round_line)]
    points <- layers[sapply(layers, is_point)]
    
    layers <- c(others, round_lines, points)
}