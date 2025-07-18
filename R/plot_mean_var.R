#' @title Plot Mean-Variance Relationship
#' @description Create a scatter plot showing the relationship between mean expression
#' and variance for genes, colored by phylostratum, with optional highlighting and labeling of specific genes.
#' 
#' @param phyex_set A PhyloExpressionSet object containing gene expression data.
#' @param highlight_genes Optional character vector of gene IDs to highlight and label on the plot.
#' @param colour_by Character string specifying coloring scheme: "none" (default), "strata" colors by phylostratum
#' 
#' @return A ggplot2 object showing the mean-variance relationship.
#' 
#' @details
#' This function plots the mean expression versus variance for each gene,
#' with points colored by phylostratum. Optionally, specific genes can be highlighted and labeled.
#' This visualization helps identify expression patterns and heteroscedasticity in the data.
#' 
#' @examples
#' # Create mean-variance plot
#' # mv_plot <- plot_mean_var(phyex_set)
#' # Highlight and label specific genes
#' # mv_plot <- plot_mean_var(phyex_set, highlight_genes = c("GeneA", "GeneB"))
#' 
#' @import ggplot2
#' @import dplyr
#' @export
plot_mean_var <- function(phyex_set, highlight_genes = NULL, colour_by = c("none", "strata")) {
    colour_by <- match.arg(colour_by)
    df <- phyex_set@data_collapsed |>
        mutate(
            mean = rowMeans(phyex_set@counts_collapsed),
            var = rowVars(phyex_set@counts_collapsed)
        ) |>
        mutate(
            mean = mean + 1e-6,
            var = var + 1e-6
        )

    # Add highlight column
    df$highlight <- if (!is.null(highlight_genes)) {
        df$GeneID %in% highlight_genes
    } else {
        FALSE
    }

    if (colour_by == "strata") {
        p <- ggplot(df, aes(x = mean, y = var)) +
            geom_point(aes(colour = Stratum), alpha = 0.7, size = 0.2) +
            scale_color_manual(values = PS_colours(length(unique(df$Stratum))), na.value = "grey50", name = "Stratum") +
            scale_x_log10() +
            scale_y_log10() +
            theme_minimal()
    } else {
        p <- ggplot(df, aes(x = mean, y = var)) +
            geom_point(colour = "black", alpha = 0.7, size = 0.2) +
            scale_x_log10() +
            scale_y_log10() +
            theme_minimal()
    }

    # Add highlighted genes in red and label them
    if (any(df$highlight)) {
        p <- p + geom_point(
            data = df[df$highlight, ],
            aes(x = mean, y = var),
            colour = "red",
            size = 1,
            alpha = 0.9,
            inherit.aes = FALSE
        ) +
        ggrepel::geom_text_repel(
            data = df[df$highlight, ],
            aes(x = mean, y = var, label = GeneID),
            colour = "red",
            size = 3,
            inherit.aes = FALSE,
            max.overlaps = 10,
            box.padding = 0.5,
            segment.color = "red",
            bg.color = "white",
            bg.r = 0.15
        )
    }

    p
}