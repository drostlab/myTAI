#' @title Comparing expression levels distributions across developmental stages
#' @description \code{plot_distribution_expression} generates plots that help to compare the distribution
#' of expression levels through various developmental stages or cell types, highlighting each stage with
#' distinct colors. By default, a log transformation is applied to the expression values.
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet).
#' @param show_identities Logical, whether to show identity-specific distributions.
#' @param show_strata Logical, whether to show stratum-specific distributions.
#' @param log_transform Logical, whether to apply log transformation to expression values (default: TRUE).
#' @param seed Seed for reproducible color selection.
#' @section Recommendation: 
#' Apply a square root transformation to enhance the visualization of differences
#' in the distributions: \code{plot_distribution_expression(transform_counts(phyex_set, sqrt))}
#' @author Filipa Martins Costa
#' @import ggplot2 tibble patchwork
#' @importFrom tidyr pivot_longer
#' @importFrom ggridges geom_density_ridges
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @export

plot_distribution_expression <- function(phyex_set,
                                         show_identities = TRUE,
                                         show_strata = FALSE,
                                         log_transform = TRUE,
                                         seed = 123) {
    # Create data frame with gene info and expression data
    expr_df <- data.frame(
        GeneID = phyex_set@gene_ids,
        Stratum = phyex_set@strata,
        phyex_set@expression_collapsed
    )
    
    df <- tidyr::pivot_longer(
        expr_df,
        cols = -c(Stratum, GeneID),
        names_to = "Identity",
        values_to = "Expression"
    )
    
    # Apply log transformation if requested
    if (log_transform) {
        df$Expression <- log1p(df$Expression)  # log1p for log(1+x) to handle zeros
    }

    qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
    col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    set.seed(seed)
    colors <- sample(col_vector, ncol(phyex_set@expression_collapsed))

    p <- ggplot(
        df,
        aes(x = Expression)
    ) +
        geom_density(alpha = 0.7, color = "black", fill = "darkgray") +
        labs(
            x = "Expression",
            y = "Density"
        ) +
        theme_minimal()

    if (show_identities) {
        p_cond <- ggplot(
            df,
            aes(x = Expression, y = factor(Identity, levels = unique(Identity)), fill = factor(Identity, levels = unique(Identity)))
        ) +
            ggridges::geom_density_ridges(alpha = 0.7, color = "black", scale = 1) +
            labs(
                x = "Expression",
                y = "Density",
                fill = phyex_set@identities_label
            ) +
            theme_minimal() +
            scale_fill_manual(values = colors)

        p <- p + p_cond
    }
    if (show_strata) {
        p_strata <- ggplot(
            df,
            aes(x = Expression, y = Stratum, fill = Stratum)
        ) +
            ggridges::geom_density_ridges(alpha = 0.7, color = "black", scale = 1) +
            labs(
                x = "Expression",
                y = "Stratum"
            ) +
            theme_minimal() +
            scale_fill_manual(values = PS_colours(phyex_set@num_strata))

        p <- p + p_strata
    }
    p + plot_layout(nrow = 1)
}
