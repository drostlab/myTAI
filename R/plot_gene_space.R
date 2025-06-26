
plot_gene_space <- function(phyex_set, top_p=0.2) {
    e <- phyex_set@counts_collapsed |>
        log1p() |>
        filter_dyn_expr(thr=1-top_p) |>
        to_std_expr()
        
        
    N = nrow(e)
    S = ncol(e)
    ideal_genes <- rbind(early_gene(S), mid_gene(S), 
                         late_gene(S), rev_mid_gene(S))
    rownames(ideal_genes) <- c("early", "mid", "late", "rev_mid")
    e_aug <- rbind(e, ideal_genes)
    
    pca <- prcomp(e_aug, scale. = FALSE)
    coords <- pca$x[, 1:2]
    
    df <- data.frame(PC1 = coords[,1], PC2 = coords[,2])
    df$label <- rownames(df)
    
    df$angle <- atan2(df$PC2, df$PC1)
    
    df$angle <- mod_pi(df$angle - df["rev_mid", "angle"] + pi)
    
    
    if (df["early", "angle"] < df["late", "angle"])
        df$angle <- mod_pi(-df$angle)
    
    ggplot(df[1:N, ], aes(PC1, PC2, color = angle)) +
        geom_point(alpha = 0.7) +
        scale_color_viridis_c() +
        geom_label(data = df[df$label %in% rownames(ideal_genes), ],
                   aes(label=label), vjust=-1, color="red") +
        coord_fixed() +
        theme_minimal() +
        ggtitle("Gene space PCA") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey")
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
                       color = colorRampPalette(c("#0055A4", "#FFFFFF", "#EF4135"))(99),
                       main="Gene expression")
}

early_gene <- function(S) {
    c(rep(-1, length.out=S%/%4),
      seq(from=-1, to=1, length.out=S-2*(S%/%4)),
      rep(1, length.out=S%/%4))
}

mid_gene <- function(S) {
    c(rep(-1, length.out=S%/%6),
      seq(from=-1, to=1, length.out=S%/%4),
      rep(1, length.out=S-2*(S%/%6)-2*(S%/%4)),
      seq(from=1, to=-1, length.out=S%/%4),
      rep(-1, length.out=S%/%6))
}

late_gene <- function(S) {
    -early_gene(S)
}

rev_mid_gene <- function(S) {
    -mid_gene(S)
}

mod_pi <- function(x) {
    (x + pi) %% (2*pi) - pi
}

# pick top variance genes from quantile
filter_dyn_expr <- function(e, thr=0.9) {
    var_genes <- apply(e, 1, var)
    cutoff <- quantile(var_genes, thr)
    e[var_genes > cutoff, ]
}

# standardise expression of genes
to_std_expr <- function(e) {
    row_sd <- apply(e, 1, sd, na.rm = TRUE)
    valid <- row_sd > 0 & is.finite(row_sd)
    e[valid, ] <- t(scale(t(e[valid, , drop = FALSE]), center = TRUE, scale = TRUE))
    e[!valid, ] <- 0
    e
}

get_angles <- function(e) {
    N = nrow(e)
    S = ncol(e)
    ideal_genes <- rbind(early_gene(S), mid_gene(S), 
                         late_gene(S), rev_mid_gene(S))
    rownames(ideal_genes) <- c("early", "mid", "late", "rev_mid")
    e_aug <- rbind(e, ideal_genes)
    
    pca <- prcomp(e_aug, scale. = FALSE)
    coords <- pca$x[, 1:2]
    
    df <- data.frame(PC1 = coords[,1], PC2 = coords[,2])
    
    df$angle <- atan2(df$PC2, df$PC1)
    
    df$angle <- mod_pi(df$angle - df["rev_mid", "angle"] + pi)
    
    if (df["early", "angle"] < df["late", "angle"])
        df$angle <- mod_pi(-df$angle)
    
    df[1:N, "angle"]
}