

#' @import ggplot2 dplyr tidyr
#' @importFrom ggrepel geom_text_repel
#' @export
plot_gene_profiles <- function(phyex_set,
                               genes = NULL,
                               show_set_mean=FALSE,
                               show_reps = FALSE,
                               transformation=c("log", "std_log", "none"),
                               colour_by=c("stage", "strata", "manual"),
                               colours = NULL,
                               show_bg = TRUE,
                               bg_top_p = 0.02,
                               show_labels=TRUE,
                               show_legend=TRUE) {
    
    transformation <- match.arg(transformation)
    colour_by <- match.arg(colour_by)
    
    
    counts <- phyex_set@counts
    counts <- switch(transformation,
                     log = log1p(counts),
                     std_log = to_std_expr(log1p(counts)),
                     none = counts)
    
    all_genes <- phyex_set@gene_ids
    
    if (show_bg)
        bg_genes <- filter_dyn_expr(counts, thr = 1 - bg_top_p) |> rownames()
    else
        bg_genes <- c()
    
    # show the gene if it is highlighted and the backgrounds if specified
    show_gene <- all_genes %in% genes | all_genes %in% bg_genes
    counts <- counts[show_gene, , drop=FALSE]
    
    df_long <- reshape2::melt(counts)
    colnames(df_long) <- c("GeneID", "Sample", "Expression")
    
    df_long <- df_long |>
        left_join(data.frame(Sample = phyex_set@sample_names,
                             Condition = phyex_set@groups),
                  by = "Sample") |>
        group_by(GeneID, Condition) |>
        summarise(min = min(Expression),
                  max = max(Expression),
                  Expression = mean(Expression),
                  .groups = "drop") |>
        mutate(Highlight = factor(GeneID %in% genes, levels = c(FALSE, TRUE))) |>
        left_join(data.frame(GeneID = phyex_set@gene_ids,
                             Stratum = phyex_set@strata,
                             Angle = -get_angles(phyex_set@counts |> log1p() |> to_std_expr())),
                  by = "GeneID") |> 
        mutate(
            ColourVar = switch(colour_by,
                               stage = Angle,
                               strata = Stratum,
                               manual = GeneID))
    
    
    p <- ggplot(df_long, aes(x = Condition,
                             y = Expression,
                             group = GeneID,
                             colour = ColourVar))  +
        geom_line(aes(alpha = Highlight, linewidth = Highlight)) +
        scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
        scale_linewidth_manual(values = c("TRUE" = 0.8, "FALSE" = 0.2), guide = "none")
        labs(x = phyex_set@conditions_label, y = "Expression") +
        theme_minimal()
    
    # show ribbon of replicates for each highlighted gene
    if (show_reps)
        p <- p + geom_ribbon(data = df_long |> filter(Highlight == 'TRUE'),
                        aes(x = Condition, ymin = min, ymax = max, 
                            fill = ColourVar, 
                            group = GeneID),
                        alpha=0.25, inherit.aes = FALSE)
    
    
    if (show_set_mean && length(genes) > 0) {
        df_mean <- df_long |> filter(Highlight == 'TRUE') |> group_by(Condition) |>
            summarise(Expression = mean(Expression), .groups="drop")
        p <- p + geom_line(data = df_mean, aes(x = Condition, y = Expression),
                           inherit.aes = FALSE,
                           colour = "red",
                           group=1,
                           linewidth = 1.2)
    }
    
    if (show_labels && length(genes) > 0) {
        df_labels <- df_long |> filter(Highlight == 'TRUE') |> group_by(GeneID) |>
            slice_max(Expression, n=1, with_ties = FALSE) 
        
        p <- p + ggrepel::geom_text_repel(data = df_labels,
                                          aes(x = Condition, y = Expression, label = GeneID),
                                          inherit.aes = FALSE,
                                          size = 3, max.overlaps=20, 
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
        gene_levels <- intersect(genes, df_long$GeneID)
        if (is.null(colours)) {
            n <- length(gene_levels)
            base_palette <- RColorBrewer::brewer.pal(min(n, 9), "Set1")
            colours <- setNames(colorRampPalette(base_palette)(n), gene_levels)
        }
        p <- p + scale_colour_manual(values = colours, name = "GeneID")
        if (show_reps)
            p <- p + scale_fill_manual(values = colours, name = "GeneID", guide = "none")
    }
    
    if (!show_legend)
        p <- p + theme(legend.position = "none")
    
    p
    
} 
