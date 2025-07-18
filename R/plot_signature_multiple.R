


#' @title Plot Multiple Transcriptomic Signatures
#' @description Create a plot comparing multiple transcriptomic signatures on the same axes,
#' with options for statistical testing and transformations.
#' 
#' @param phyex_sets A vector of PhyloExpressionSet objects to compare
#' @param legend_title Title for the legend (default: "Phylo Expression Set")
#' @param show_p_val Logical indicating whether to show p-values (default: FALSE)
#' @param conservation_test Function to use for conservation testing (default: flatline_test)
#' @param transformation Optional transformation function to apply to all datasets (default: NULL)
#' @param colours Optional vector of colors for each dataset (default: NULL)
#' @param ... Additional arguments passed to plot_signature
#' 
#' @return A ggplot2 object showing multiple transcriptomic signatures
#' 
#' @details
#' This function allows comparison of multiple transcriptomic signatures by overlaying
#' them on the same plot. Each signature is colored differently and can be tested
#' for conservation patterns. If a transformation is provided, it's applied to all
#' datasets before plotting.
#' 
#' @examples
#' # Compare multiple datasets
#' # phyex_list <- c(phyex_set1, phyex_set2)
#' # p <- plot_signature_multiple(phyex_list, legend_title = "Dataset")
#' 
#' # With transformation
#' # p2 <- plot_signature_multiple(phyex_list, transformation = log1p)
#' 
#'
#' @import ggplot2 purrr
#' @export
plot_signature_multiple <- function(phyex_sets,
                                    legend_title="Phylo Expression Set",
                                    show_p_val=FALSE,
                                    conservation_test=flatline_test,
                                    transformation = NULL,
                                    colours=NULL,
                                    ...) {
    
    if (!is.null(transformation))
        phyex_sets <- phyex_sets |>
            map(tf, FUN=transformation)
    layers <- phyex_sets |>
        map(plot_signature, show_p_val=F, ...) |>
        map(~ .x$layers) |>
        unlist(..., recursive = TRUE) |>
        .sort_layers()
    
    p <- ggplot() +
        layers +
        labs(
            x=phyex_sets[[1]]@conditions_label,
            y=phyex_sets[[1]]@index_full_name,
            colour=legend_title,
            fill=legend_title
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    labels <- phyex_sets |> map_chr(\(s) s@name)
    
    # add p values to legend
    if (show_p_val) {
        test_res <- phyex_sets |> map(conservation_test, plot_result=FALSE)
        p_vals <- test_res |> map(\(t) t@p_value) |>
            map(\(pval) signif(pval, digits=3)) |> as.character()
        p_label <- test_res[[1]]@p_label
        labels <- map2(labels, p_vals, \(l, pval) paste0(l, "<br>**", p_label, "**: ", pval, "<br>"))
        p <- p + theme(legend.text=ggtext::element_markdown())
    }
    
    if (!is.null(colours))
        p <- p + scale_color_manual(labels=labels, values=colours) + scale_fill_manual(labels=labels, values=colours)
    else
        p <- p + scale_color_hue(labels=labels) + scale_fill_hue(labels=labels)
        
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