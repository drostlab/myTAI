

#' @import ggplot2 dplyr tidyr
#' @importFrom ggrepel geom_text_repel
plot_gene_profiles <- function(phyex_set,
                               genes = NULL,
                               show_bg = FALSE,
                               show_set_mean=FALSE,
                               show_reps = FALSE,
                               transformation=c("log", "std_log", "none"),
                               colour_by=c("stage", "strata", "manual"),
                               colours = NULL,
                               top_p = 1.0,
                               max_bg_genes = 500,
                               show_labels=TRUE,
                               show_legend=TRUE) {
    
    transformation <- match.arg(transformation)
    colour_by <- match.arg(colour_by)
    
    
    counts <- phyex_set@counts
    counts <- switch(transformation,
                     log = log1p(counts),
                     std_log = to_std_expr(log1p(counts)),
                     none = counts)
    
    if (is.null(genes)) {
        filtered <- filter_dyn_expr(counts, thr = 1 - top_p)
        genes <- rownames(filtered)
    }
    
    all_genes <- phyex_set@gene_ids
    is_highlighted <- all_genes %in% genes
    # show the gene if it is highlighted and the backgrounds if specified
    show_gene <- is_highlighted | show_bg
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
                             Stratum = phyex_set@stratas,
                             Angle = -get_angles(phyex_set@counts |> log1p() |> to_std_expr())),
                  by = "GeneID") |>
        mutate(Colour = switch(colour_by,
                               stage = Angle,
                               strata = Stratum,
                               manual = GeneID))
    
    # show only a subset of the bg genes
    bg_genes <- df_long |> filter(Highlight == "FALSE") |> distinct(GeneID) |> pull(GeneID)
    set.seed(123)
    sampled_bg_genes <- sample(bg_genes, size = min(max_bg_genes, length(bg_genes)))
    
    p <- ggplot(df_long, aes(x = Condition,
                             y = Expression,
                             group = GeneID,
                             colour = Colour,
                             alpha = Highlight,
                             linewidth = Highlight)) +
        geom_line(data = df_long |> filter(GeneID %in% sampled_bg_genes)) +
        geom_line(data = df_long |> filter(Highlight == "TRUE")) +
        scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
        scale_linewidth_manual(values = c("TRUE" = 0.8, "FALSE" = 0.2), guide = "none")
        labs(x = phyex_set@conditions_label, y = "Expression") +
        theme_minimal()
    
    # show ribbon of replicates for each highlighted gene
    if (show_reps)
        p <- p + geom_ribbon(data = df_long |> filter(Highlight == 'TRUE'),
                        aes(x = Condition, ymin = min, ymax = max, fill = Colour, group = GeneID),
                        alpha=0.4, inherit.aes = FALSE)
    
    
    if (show_set_mean && length(genes) > 0) {
        df_mean <- df_long |> filter(Highlight == 'TRUE') |> group_by(Condition) |>
            summarise(Expression = mean(Expression), .groups="drop")
        p <- p + geom_line(data = df_mean, aes(x = Condition, y = Expression),
                           inherit.aes = FALSE,
                           colour = "red",
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
        p <- p + scale_colour_viridis_c(name = "Angle", option = "C")
        if (show_reps)
            p <- p + scale_fill_viridis_c(name = "Angle", option = "C", guide = "none")
    } else if (colour_by == "strata") {
        p <- p + scale_colour_viridis_d(name = "Stratum")
        if (show_reps)
            p <- p + scale_fill_viridis_d(name = "Stratum", guide = "none")
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



















#' @import ggplot2
plot_gene_profilesa <- function(phyex_set,
                               show_CI=T,
                               show_mean=TRUE,
                               colours=NULL) {
    
    if (is.null(colours))
        colours = rep("black", phyex_set@num_genes)
    
    df <- phyex_set@data_collapsed |>
        tidyr::pivot_longer(-c(Stratum, GeneID), names_to="Sample", values_to="Expression") |>
        left_join(data.frame(Condition=phyex_set@groups, Sample=phyex_set@sample_names), by="Sample") |>
        group_by(Stratum, GeneID, Condition) |>
        summarise(min=min(Expression), max=max(Expression), Expression=mean(Expression), .groups="drop")
    
    
    
    
    
    p <- ggplot(df, aes(x=factor(Condition, unique(Condition)),
                        group=GeneID,
                        colour=GeneID)) +
        geom_line(aes(y=Expression), alpha=0.5) +
        scale_color_manual(values=colours) +
        labs(x=phyex_set@conditions_label,
             y="Expression") +
        scale_x_discrete(labels = ~ stringr::str_wrap(., 10)) +
        theme_minimal()
    if (show_CI)
        p <- p + 
            geom_ribbon(aes(ymin=min, ymax=max), fill="gray", alpha=0.05)
    
    if (show_mean) {
        df_mean <- df |> group_by(Condition) |> summarise(Expression=mean(Expression), .groups="drop")
        p <- p + 
            geom_line(data = df_mean,
                      aes(y=Expression), colour="red", size=2, group=0)
    }
    
    
    return(p)
}