#' @title Plot Individual Gene Expression Profiles
#' @description Create a plot showing expression profiles for individual genes across
#' developmental stages or cell types, with various visualization options.
#' 
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)
#' @param genes Character vector of gene IDs to plot. If NULL, top expressing genes are selected
#' @param show_set_mean Logical indicating whether to show the mean expression across all genes (default: FALSE)
#' @param show_reps Logical indicating whether to show individual replicates (bulk) or cells (single-cell) (default: FALSE)
#' @param transformation Character string specifying expression transformation:
#' "log" (log1p), "std_log" (standardized log1p), or "none" (default: "log")
#' @param colour_by Character string specifying coloring scheme:
#' "strata" (by phylostratum), "stage" (by developmental stage/cell type), or "manual" (default: "strata")
#' @param colours Optional vector of colors for manual coloring (default: NULL)
#' @param max_genes Maximum number of genes to plot when genes=NULL (default: 100)
#' @param show_labels Logical indicating whether to show gene labels (default: TRUE)
#' @param show_legend Logical indicating whether to show legend (default: TRUE)
#' @param facet_by_strata Logical indicating whether to facet by phylostratum (default: FALSE)
#' 
#' @return A ggplot2 object showing gene expression profiles
#' 
#' @details
#' This function creates detailed visualizations of individual gene expression patterns
#' across development (bulk data) or cell types (single-cell data). Genes can be colored by phylostratum or developmental stage,
#' and various transformations can be applied to highlight different aspects of the data.
#' 
#' @examples
#' # Plot specific genes for bulk data
#' # p1 <- plot_gene_profiles(bulk_phyex_set, genes = c("gene1", "gene2", "gene3"))
#' 
#' # Plot for single-cell data with faceting by strata
#' # p2 <- plot_gene_profiles(sc_phyex_set, facet_by_strata = TRUE)
#' # p2 <- plot_gene_profiles(phyex_set, facet_by_strata = TRUE)
#' 
#' @author Filipa Martins Costa, Stefan Manolache, Hajk-Georg Drost
#' 
#' @import ggplot2 dplyr tidyr
#' @importFrom ggrepel geom_text_repel
#' @export
plot_gene_profiles <- function(phyex_set,
                               genes = NULL,
                               show_set_mean=FALSE,
                               show_reps = FALSE,
                               transformation=c("log", "std_log", "none"),
                               colour_by=c("strata", "stage", "manual"),
                               colours = NULL,
                               max_genes = 100,
                               show_labels=TRUE,
                               show_legend=TRUE,
                               facet_by_strata = FALSE) {
    
    transformation <- match.arg(transformation)
    colour_by <- match.arg(colour_by)
    
    
    counts <- phyex_set@counts
    counts <- switch(transformation,
                     log = log1p(counts),
                     std_log = to_std_expr(log1p(counts)),
                     none = counts)
    
    all_genes <- phyex_set@gene_ids
    
    # Select genes to plot
    if (is.null(genes)) {
        # If no specific genes provided, select top expressing genes
        gene_means <- rowMeans(counts)
        genes_to_plot <- names(sort(gene_means, decreasing = TRUE))[1:min(max_genes, length(gene_means))]
    } else {
        # Use only the provided genes
        genes_to_plot <- genes
    }
    
    # Filter to genes to plot
    show_gene <- all_genes %in% genes_to_plot
    counts <- counts[show_gene, , drop=FALSE]
    
    df_long <- reshape2::melt(counts)
    colnames(df_long) <- c("GeneID", "Sample", "Expression")
    
    df_long <- df_long |>
        left_join(data.frame(Sample = phyex_set@sample_names,
                             Identity = phyex_set@groups),
                  by = "Sample") |>
        group_by(GeneID, Identity) |>
        summarise(min = min(Expression),
                  max = max(Expression),
                  Expression = mean(Expression),
                  .groups = "drop") |>
        left_join(data.frame(GeneID = phyex_set@gene_ids,
                             Stratum = phyex_set@strata,
                             Angle = -get_angles(phyex_set@expression |> log1p() |> to_std_expr())),
                  by = "GeneID") |> 
        mutate(
            ColourVar = switch(colour_by,
                               stage = Angle,
                               strata = Stratum,
                               manual = GeneID))
    
    
    p <- ggplot(df_long, aes(x = Identity,
                             y = Expression,
                             group = GeneID,
                             colour = ColourVar))  +
        geom_line() +
        labs(x = phyex_set@identities_label, y = "Expression") +
        theme_minimal()
    
    # show ribbon of replicates
    if (show_reps)
        p <- p + geom_ribbon(aes(x = Identity, ymin = min, ymax = max, 
                                fill = ColourVar, 
                                group = GeneID),
                            alpha=0.25, inherit.aes = FALSE)
    
    
    if (show_set_mean && length(genes_to_plot) > 0) {
        df_mean <- df_long |> group_by(Identity) |>
            summarise(Expression = mean(Expression), .groups="drop")
        p <- p + geom_line(data = df_mean, aes(x = Identity, y = Expression),
                           inherit.aes = FALSE,
                           colour = "red",
                           group=1,
                           linewidth = 1.2)
    }
    
    if (show_labels && length(genes_to_plot) > 0) {
        if (facet_by_strata) {
            # When faceting, choose top labels within each stratum
            df_labels <- df_long |> group_by(GeneID) |>
                slice_max(Expression, n=1, with_ties = FALSE) |>
                ungroup() |>
                group_by(Stratum) |>
                slice_max(Expression, n=5, with_ties = FALSE) |>
                ungroup()
        } else {
            # Choose top labels globally
            df_labels <- df_long |> group_by(GeneID) |>
                slice_max(Expression, n=1, with_ties = FALSE) |>
                ungroup() |>
                slice_max(Expression, n=10, with_ties = FALSE)
        }
        
        p <- p + ggrepel::geom_text_repel(data = df_labels,
                                          aes(x = Identity, y = Expression, label = GeneID),
                                          inherit.aes = FALSE,
                                          size = 1.8, max.overlaps=20, 
                                          box.padding = 0.5, point.padding = 0.5,
                                          segment.size = 0.2, segment.colour = "grey50")
    }
    
    
    # handle colouring
    if (colour_by == "stage") {
        p <- p + scale_colour_viridis_c(name = "Angle")
        if (show_reps)
            p <- p + scale_fill_viridis_c(name = "Angle", guide = "none")
    } else if (colour_by == "strata") {
        levels <- levels(phyex_set@strata)
        p <- p + scale_color_manual(name = "Strata", values = PS_colours(phyex_set@num_strata), na.value = "grey50", limits=levels, drop=FALSE)
        if (show_reps)
            p <- p + scale_fill_manual(name = "Strata", values = PS_colours(phyex_set@num_strata), na.value = "grey50", limits=levels, drop=FALSE)
    } else if (colour_by == "manual") {
        gene_levels <- intersect(genes_to_plot, df_long$GeneID)
        if (is.null(colours)) {
            n <- length(gene_levels)
            base_palette <- RColorBrewer::brewer.pal(min(n, 9), "Set1")
            colours <- stats::setNames(grDevices::colorRampPalette(base_palette)(n), gene_levels)
        }
        p <- p + scale_colour_manual(values = colours, name = "GeneID")
        if (show_reps)
            p <- p + scale_fill_manual(values = colours, name = "GeneID", guide = "none")
    }
    
    # Add faceting if requested
    if (facet_by_strata) {
        p <- p + 
            facet_wrap(~ Stratum, scales = "fixed") +
            labs(x = paste(phyex_set@identities_label, "by Stratum"))
        
        if (colour_by == "strata") {
            p <- p + guides(colour = "none", fill = "none")
        }
    }
    
    # Add slanted x-axis labels
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    if (!show_legend)
        p <- p + theme(legend.position = "none")
    
    p
    
}
