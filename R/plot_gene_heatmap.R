#' @title Plot Gene Expression Heatmap
#' @description Create a heatmap showing gene expression patterns across conditions
#' with optional dendrograms and gene age annotation.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param genes Character vector of specific genes to plot. If NULL, uses top dynamic genes
#' @param top_p Proportion of most dynamic genes to include (default: 0.2)
#' @param std Logical indicating whether to use standardized expression values (default: TRUE)
#' @param reps Logical indicating whether to show individual replicates or collapsed conditions (default: FALSE)
#' @param cluster_rows Logical indicating whether to cluster genes/rows (default: FALSE)
#' @param cluster_cols Logical indicating whether to cluster conditions/columns (default: FALSE)
#' @param show_gene_age Logical indicating whether to show gene age as row annotation (default: TRUE)
#' @param show_gene_ids Logical indicating whether to show gene names (default: FALSE)
#' @param ... Additional arguments passed to pheatmap::pheatmap
#' 
#' @return A ggplot object (converted from pheatmap) showing the gene expression heatmap
#' 
#' @details
#' This function creates a comprehensive heatmap visualization of gene expression patterns.
#' By default, genes are ordered by their expression angle (developmental trajectory).
#' The function supports clustering of both genes and conditions, and can optionally
#' display gene age (phylostratum) as a colored annotation bar.
#' 
#' The gene age annotation uses the PS_colours function to create a consistent
#' color scheme across different myTAI visualizations.
#' 
#' @examples
#' # Basic heatmap with gene age annotation
#' # p1 <- plot_gene_heatmap(phyex_set, show_gene_age = TRUE)
#' 
#' # Clustered heatmap with standardized values
#' # p2 <- plot_gene_heatmap(phyex_set, std = TRUE, cluster_rows = TRUE, cluster_cols = TRUE)
#' 
#' # Heatmap of specific genes with replicates
#' # p3 <- plot_gene_heatmap(phyex_set, genes = c("gene1", "gene2"), reps = TRUE)
#' 
#' @export
plot_gene_heatmap <- function(phyex_set, 
                              genes = NULL,
                              top_p = 0.2, 
                              std = TRUE, 
                              reps = FALSE,
                              cluster_rows = FALSE,
                              cluster_cols = FALSE,
                              show_gene_age = TRUE,
                              show_gene_ids = FALSE,
                              ...) {

    
    # Select expression data
    if (reps) {
        e <- phyex_set@counts
    } else {
        e <- phyex_set@counts_collapsed
    }
    
    # Apply log transformation
    e <- e |> log1p()
    
    # Filter genes if specific genes are provided
    if (!is.null(genes) && length(genes) > 0) {
        e <- e[rownames(e) %in% genes, , drop = FALSE]
    } else {
        # Filter for dynamic genes
        e <- e |> filter_dyn_expr(thr = 1 - top_p)
    }
    
    # Calculate standardized expression for ordering
    se <- e |> to_std_expr()
    
    # Use standardized values if requested
    if (std) {
        e <- se
    }
    
    # Order genes by expression angle if not clustering
    if (!cluster_rows) {
        gene_order <- order(get_angles(se))
        e <- e[gene_order, ]
    }
    
    
    color_palette <- grDevices::colorRampPalette(c("#0055A4", "#FFFFFF", "#EF4135"))(99)
        
    
    # Prepare annotations
    annotation_row <- NULL
    annotation_colors <- NULL
    
    if (show_gene_age) {
        # Create gene age annotation
        gene_names <- rownames(e)
        strata_map <- stats::setNames(phyex_set@strata, phyex_set@gene_ids)
        gene_strata <- strata_map[gene_names]
        
        # Remove genes not found in the mapping
        gene_strata <- gene_strata[!is.na(gene_strata)]
        
        # Create annotation data frame
        annotation_row <- data.frame(
            Phylostratum = gene_strata,
            row.names = names(gene_strata)
        )
        
        # Create color mapping for phylostrata using all strata levels
        all_strata_levels <- levels(phyex_set@strata)
        ps_colors <- PS_colours(phyex_set@num_strata)
        names(ps_colors) <- all_strata_levels
        
        annotation_colors <- list(
            Phylostratum = ps_colors
        )
    }
    
    p <- pheatmap::pheatmap(
        e,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_rownames = show_gene_ids,
        show_colnames = TRUE,
        color = color_palette,
        annotation_row = annotation_row,
        annotation_colors = annotation_colors,
        fontsize_row = 5,
        silent = TRUE,
        ...
    ) |>
        ggplotify::as.ggplot()
    
    # Convert to ggplot
    return(p)
}