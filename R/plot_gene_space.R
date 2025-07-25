
#' @title Plot Gene Space Using PCA
#' @description Create a PCA plot showing genes in expression space with ideal expression
#' patterns overlaid as reference points.
#' 
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)
#' @param top_p Proportion of most dynamic genes to include when genes=NULL (default: 0.2)
#' @param genes Character vector of specific genes to plot. If NULL, uses top dynamic genes
#' @param colour_by Character string specifying coloring scheme:
#' "identity" (by peak expression stage/identity) or "strata" (by phylostratum) (default: "identity")
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
#' # Plot gene space colored by identity
#' # p1 <- plot_gene_space(phyex_set, colour_by = "identity")
#' 
#' # Plot specific genes colored by strata
#' # p2 <- plot_gene_space(phyex_set, genes = c("gene1", "gene2"), colour_by = "strata")
#' 
#' @export
plot_gene_space <- function(phyex_set, 
                            top_p=0.2,
                            genes=NULL,
                            colour_by=c("identity", "strata")) {
    colour_by <- match.arg(colour_by)
    
    e <- phyex_set@expression_collapsed |>
        log1p()
    
    if (!is.null(genes) && length(genes) > 0) {
        e <- e[rownames(e) %in% genes, , drop=FALSE]
    } else {
        e <- e |> genes_filter_dynamic(thr = 1 - top_p)
    }
    
    e <- .to_std_expr(e)
        
        
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
    
    
    
    p <- ggplot(df[1:N, ], aes(PC1, PC2, color = if (colour_by == "identity") angle else strata)) +
        geom_point(aes(alpha=highlight)) +
        scale_alpha_manual(values = c(`TRUE` = 0.7, `FALSE` = 0.1), guide = "none")
    if (colour_by == "identity") {
        p <- p + 
            scale_color_viridis_c(name = phyex_set@identities_label)
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