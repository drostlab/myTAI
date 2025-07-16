
#' @title Plot Gene Space Using PCA
#' @description Create a PCA plot showing genes in expression space with ideal expression
#' patterns overlaid as reference points.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param top_p Proportion of most dynamic genes to include when genes=NULL (default: 0.2)
#' @param genes Character vector of specific genes to plot. If NULL, uses top dynamic genes
#' @param colour_by Character string specifying coloring scheme:
#' "stage" (by peak expression stage) or "strata" (by phylostratum) (default: "stage")
#' 
#' @return A ggplot2 object showing the gene space PCA plot
#' 
#' @details
#' This function creates a PCA visualization of genes in expression space, with ideal
#' expression patterns (early, mid, late, reverse mid) overlaid as reference points.
#' The analysis uses log-transformed and standardized expression values. Genes are
#' colored either by their phylostratum or by their peak expression stage.
#' 
#' @examples
#' # Plot gene space colored by stage
#' # p1 <- plot_gene_space(phyex_set, colour_by = "stage")
#' 
#' # Plot specific genes colored by strata
#' # p2 <- plot_gene_space(phyex_set, genes = c("gene1", "gene2"), colour_by = "strata")
#' 
#' @export
plot_gene_space <- function(phyex_set, 
                            top_p=0.2,
                            genes=NULL,
                            colour_by=c("stage", "strata")) {
    colour_by <- match.arg(colour_by)
    
    e <- phyex_set@counts_collapsed |>
        log1p()
    
    if (!is.null(genes) && length(genes) > 0) {
        e <- e[rownames(e) %in% genes, , drop=FALSE]
    } else {
        e <- e |> filter_dyn_expr(thr = 1 - top_p)
    }
    
    e <- to_std_expr(e)
        
        
    N = nrow(e)
    S = ncol(e)
    ideal_genes <- rbind(early_gene(S), mid_gene(S), 
                         late_gene(S), rev_mid_gene(S))
    rownames(ideal_genes) <- c("early", "mid", "late", "rev_mid")
    e_aug <- rbind(e, ideal_genes)
    
    pca <- stats::prcomp(e_aug, scale. = FALSE)
    coords <- pca$x[, 1:2]
    
    df <- data.frame(PC1 = coords[,1], PC2 = coords[,2])
    df$label <- rownames(df)
    
    df$angle <- atan2(df$PC2, df$PC1)
    
    df$angle <- mod_pi(df$angle - df["rev_mid", "angle"] + pi)
    
    if (df["early", "angle"] < df["late", "angle"])
        df$angle <- mod_pi(-df$angle)
    
    df$highlight <- if (is.null(genes)) TRUE else df$label %in% genes
    
    strata_map <- stats::setNames(phyex_set@strata, phyex_set@gene_ids)
    df$strata <- strata_map[df$label]
    
    
    
    p <- ggplot(df[1:N, ], aes(PC1, PC2, color = if (colour_by == "stage") angle else strata)) +
        geom_point(aes(alpha=highlight)) +
        scale_alpha_manual(values = c(`TRUE` = 0.7, `FALSE` = 0.1), guide = "none")
    if (colour_by == "stage") {
        p <- p + 
            scale_color_viridis_c(name = "Stage")
    } else {
        p <- p + 
            scale_color_manual(values = PS_colours(length(unique(df$strata))), na.value = "grey50", name = "Stratum")
    }
        
    p <- p +
        geom_label(data = df[df$label %in% rownames(ideal_genes), ],
                   aes(label = label), vjust = -1, color = "red") +
        coord_fixed() +
        theme_minimal() +
        ggtitle("Gene space PCA") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey")     
        
    p
}

plot_gene_profiles_sorted <- function(phyex_set, top_p=0.05, std=FALSE, reps=FALSE, highlighted_gene=NULL) {
    if (reps)
        e <- phyex_set@counts
    else
        e <- phyex_set@counts_collapsed
    e <- e |>
        log1p() |>
        filter_dyn_expr(thr=1-top_p)
    
    se <- e |> to_std_expr()
    
    if (std)
        e <- se
    df_long <- reshape2::melt(e)
    colnames(df_long) <- c("Gene", "Sample", "Expression")
    
    df <- data.frame(Gene = rownames(e), Angle = -get_angles(se))
    df_long <- merge(df_long, df, by = "Gene")
    
    p <- ggplot(df_long, aes(Sample, Expression, group = Gene, colour = Angle)) +
        geom_line(alpha = 1 - sqrt(top_p)) +
        labs(x=phyex_set@conditions_label) +
        scale_color_viridis_c() +
        theme_minimal()
    
    if (!is.null(highlighted_gene) && highlighted_gene %in% df_long$Gene) {
        df_highlight <- df_long[df_long$Gene == highlighted_gene, ]
        print(df_highlight) 
        p <- p + geom_line(data = df_highlight, aes(Sample, Expression, group = Gene),
                           colour = "red", linewidth = 1.2)
    }
    
    p
}

plot_gene_heatmap <- function(phyex_set, top_p=0.2, std=FALSE, reps=FALSE) {
    if (reps)
        e <- phyex_set@counts
    else
        e <- phyex_set@counts_collapsed
    e <- e |>
        log1p() |>
        filter_dyn_expr(thr=1-top_p)
    se <- e |> to_std_expr()
    if (std)
        e <- se
    
    gene_order <- order(get_angles(se))
    
    ordered_expr <- e[gene_order, ]
    pheatmap::pheatmap(ordered_expr, cluster_rows = FALSE, cluster_cols = FALSE,
                       show_rownames = FALSE,
                       color = grDevices::colorRampPalette(c("#0055A4", "#FFFFFF", "#EF4135"))(99),
                       main="Gene expression")
}

