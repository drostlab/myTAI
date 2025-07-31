#' @title Plot Gene Expression Heatmap
#' @description Create a heatmap showing gene expression patterns across conditions
#' with optional dendrograms and gene age annotation.
#' 
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)
#' @param genes Character vector of specific gene IDs to include (default: NULL for auto-selection)
#' @param top_p Numeric value specifying the top proportion of genes to include (default: 0.2)
#' @param std Logical indicating whether to standardize expression values (default: TRUE)
#' @param show_reps Logical indicating whether to show replicates or collapsed data (default: FALSE)
#' @param cluster_rows Logical indicating whether to cluster genes/rows (default: FALSE)
#' @param cluster_cols Logical indicating whether to cluster conditions/columns (default: FALSE)
#' @param show_gene_age Logical indicating whether to show gene age annotation (default: TRUE)
#' @param show_gene_ids Logical indicating whether to show gene identifiers (default: FALSE)
#' @param ... Additional arguments passed to specific methods
#' 
#' @return A ggplot object (converted from pheatmap) showing the gene expression heatmap
#' 
#' @details
#' This function creates a comprehensive heatmap visualization of gene expression patterns.
#' By default, genes are ordered by their expression angle (developmental trajectory).
#' The function supports clustering of both genes and identities, and can optionally
#' display gene age (phylostratum) as a colored annotation bar.
#' 
#' For bulk data, the heatmap shows expression across developmental conditions.
#' For single-cell data, the heatmap shows expression across cell types.
#' 
#' The gene age annotation uses the PS_colours function to create a consistent
#' color scheme across different myTAI visualizations.
#' 
#' @examples
#' # Basic heatmap with gene age annotation
#' # p1 <- plot_gene_heatmap(bulk_phyex_set, show_gene_age = TRUE)
#' 
#' # Single-cell heatmap with subset of cells
#' # p2 <- plot_gene_heatmap(sc_phyex_set, show_reps = TRUE, max_cells_per_type = 3)
#' 
#' @export
plot_gene_heatmap <- S7::new_generic("plot_gene_heatmap", "phyex_set",
    function(phyex_set,
             genes = NULL,
             top_p = 0.2,
             std = TRUE,
             show_reps = FALSE,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             show_gene_age = TRUE,
             show_gene_ids = FALSE,
             ...) {
        S7::S7_dispatch()
    }
)

#' @title Shared Gene Heatmap Implementation
#' @description Internal helper function that contains the shared logic for creating gene heatmaps.
#' 
#' @param expression_matrix Matrix of expression values (genes x samples)
#' @param strata Factor vector of gene phylostrata
#' @param gene_ids Character vector of gene IDs
#' @param num_strata Integer number of phylostrata
#' @param genes Character vector of specific genes to plot. If NULL, uses top dynamic genes
#' @param top_p Proportion of most dynamic genes to include (default: 0.2)
#' @param std Logical indicating whether to use standardized expression values (default: TRUE)
#' @param cluster_rows Logical indicating whether to cluster genes/rows (default: FALSE)
#' @param cluster_cols Logical indicating whether to cluster identities/columns (default: FALSE)
#' @param show_gene_age Logical indicating whether to show gene age as row annotation (default: TRUE)
#' @param show_gene_ids Logical indicating whether to show gene names (default: FALSE)
#' @param ... Additional arguments passed to pheatmap::pheatmap
#' 
#' @return A ggplot object showing the gene expression heatmap
#' 
#' @keywords internal
.plot_gene_heatmap_impl <- function(expression_matrix, 
                                   strata, 
                                   gene_ids, 
                                   num_strata,
                                   genes = NULL,
                                   top_p = 0.2, 
                                   std = TRUE, 
                                   cluster_rows = FALSE,
                                   cluster_cols = FALSE,
                                   show_gene_age = TRUE,
                                   show_gene_ids = FALSE,
                                   ...) {
    
    # Apply log transformation
    e <- expression_matrix |> log1p()
    
    # Filter genes if specific genes are provided
    if (!is.null(genes) && length(genes) > 0) {
        e <- e[rownames(e) %in% genes, , drop = FALSE]
    } else {
        # Filter for dynamic genes
        e <- e |> genes_filter_dynamic(thr = 1 - top_p)
    }
    
    # Calculate standardized expression for ordering
    se <- e |> .to_std_expr()
    
    # Use standardized values if requested
    if (std) {
        e <- se
    }
    

    # If only one gene, disable clustering
    n_genes <- nrow(e)
    cluster_rows <- cluster_rows
    if (n_genes == 1) {
        cluster_rows <- FALSE
    }
    # Order genes by expression angle if not clustering
    if (!cluster_rows) {
        gene_order <- order(get_angles(se))
        e <- e[gene_order, , drop = FALSE]
    }
    
    color_palette <- grDevices::colorRampPalette(c("#0055A4", "#FFFFFF", "#EF4135"))(99)
    
    # Prepare annotations
    annotation_row <- NULL
    annotation_colors <- NULL
    
    if (show_gene_age) {
        # Create gene age annotation
        gene_names <- rownames(e)
        strata_map <- stats::setNames(strata, gene_ids)
        gene_strata <- strata_map[gene_names]
        # Remove genes not found in the mapping
        gene_strata <- gene_strata[!is.na(gene_strata)]
        # Create annotation data frame (only present genes)
        annotation_row <- data.frame(
            Phylostratum = gene_strata,
            row.names = names(gene_strata)
        )
        # Always use all strata levels for color mapping
        all_strata_levels <- levels(strata)
        ps_colors <- PS_colours(num_strata)
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
    
    return(p)
}

#' @export
S7::method(plot_gene_heatmap, BulkPhyloExpressionSet) <- function(phyex_set,
                                                                 genes = NULL,
                                                                 top_p = 0.2, 
                                                                 std = TRUE, 
                                                                 show_reps = FALSE,
                                                                 cluster_rows = FALSE,
                                                                 cluster_cols = FALSE,
                                                                 show_gene_age = TRUE,
                                                                 show_gene_ids = FALSE,
                                                                 ...) {
    
    # Select expression data
    if (show_reps) {
        expression_matrix <- phyex_set@expression
    } else {
        expression_matrix <- phyex_set@expression_collapsed
    }
    
    # Call shared implementation
    .plot_gene_heatmap_impl(
        expression_matrix = expression_matrix,
        strata = phyex_set@strata,
        gene_ids = phyex_set@gene_ids,
        num_strata = phyex_set@num_strata,
        genes = genes,
        top_p = top_p,
        std = std,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_gene_age = show_gene_age,
        show_gene_ids = show_gene_ids,
        ...
    )
}

#' @export
S7::method(plot_gene_heatmap, ScPhyloExpressionSet) <- function(phyex_set,
                                                               genes = NULL,
                                                               top_p = 0.2, 
                                                               std = TRUE, 
                                                               show_reps = FALSE,
                                                               max_cells_per_type = 5,
                                                               cluster_rows = FALSE,
                                                               cluster_cols = FALSE,
                                                               show_gene_age = TRUE,
                                                               show_gene_ids = FALSE,
                                                               ...) {
    
    # Select expression data
    if (show_reps) {
        # Use downsample_expression to get a dense matrix of sampled cells
        expression_matrix <- downsample_expression(phyex_set, downsample = max_cells_per_type)
        # Warn about subsampling
        total_cells <- length(phyex_set@groups)
        n_selected <- ncol(expression_matrix)
        if (n_selected < total_cells) {
            message(paste("Showing", n_selected, "out of", total_cells, "cells.",
                         "Use max_cells_per_type to control the number of cells per type."))
        }
    } else {
        expression_matrix <- phyex_set@expression_collapsed
    }

    # Call shared implementation
    .plot_gene_heatmap_impl(
        expression_matrix = expression_matrix,
        strata = phyex_set@strata,
        gene_ids = phyex_set@gene_ids,
        num_strata = phyex_set@num_strata,
        genes = genes,
        top_p = top_p,
        std = std,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_gene_age = show_gene_age,
        show_gene_ids = show_gene_ids,
        ...
    )
}