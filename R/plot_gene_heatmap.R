#' @title Plot Gene Expression Heatmap
#' @description Create a heatmap showing gene expression patterns across conditions
#' with optional dendrograms and gene age annotation.
#' 
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)
#' @param genes Character vector of specific gene IDs to include in the heatmap (default: NULL for auto-selection of dynamic genes)
#' @param top_p Numeric value specifying the top proportion of genes to include (default: 0.2)
#' @param std Logical indicating whether to standardize expression values (default: TRUE)
#' @param show_reps Logical indicating whether to show replicates or collapsed data (default: FALSE)
#' @param cluster_rows Logical indicating whether to cluster genes/rows (default: FALSE)
#' @param cluster_cols Logical indicating whether to cluster conditions/columns (default: FALSE)
#' @param show_gene_age Logical indicating whether to show gene age annotation (default: TRUE)
#' @param show_gene_ids Logical indicating whether to show gene identifiers (default: FALSE)
#' @param gene_annotation Data frame with custom gene annotations, rownames should match gene IDs (default: NULL)
#' @param gene_annotation_colors Named list of color vectors for custom gene annotations (default: NULL)
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
#' Custom gene annotations can be provided via the \code{gene_annotation} parameter,
#' which should be a data frame with gene IDs as rownames and annotation categories
#' as columns. Corresponding colors should be provided via \code{gene_annotation_colors}
#' as a named list where names match the annotation column names.
#' 
#' @examples
#' # Basic heatmap with gene age annotation
#' # p1 <- plot_gene_heatmap(bulk_phyex_set, show_gene_age = TRUE)
#' 
#' # Single-cell heatmap with subset of cells
#' # p2 <- plot_gene_heatmap(sc_phyex_set, show_reps = TRUE, max_cells_per_type = 3)
#' 
#' # Custom gene annotation example
#' # gene_annot <- data.frame(
#' #   Category = c("High", "Medium", "Low"),
#' #   row.names = c("Gene1", "Gene2", "Gene3")
#' # )
#' # colors <- list(Category = c("High" = "red", "Medium" = "yellow", "Low" = "blue"))
#' # p3 <- plot_gene_heatmap(phyex_set, gene_annotation = gene_annot, 
#' #                        gene_annotation_colors = colors, show_gene_age = FALSE)
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
             gene_annotation = NULL,
             gene_annotation_colors = NULL,
             ...) {
        S7::S7_dispatch()
    }
)

#' @title Shared Gene Heatmap Implementation
#' @description Internal helper function that contains the shared logic for creating gene heatmaps.
#' 
#' @param expression_matrix Matrix of expression values (genes x samples)
#' @param strata Factor vector of gene phylostrata
#' @param gene_ids Character vector of all gene IDs in the dataset (used for phylostratum mapping)
#' @param num_strata Integer number of phylostrata
#' @param genes Character vector of specific genes to plot. If NULL, uses top dynamic genes
#' @param top_p Proportion of most dynamic genes to include (default: 0.2)
#' @param std Logical indicating whether to use standardized expression values (default: TRUE)
#' @param cluster_rows Logical indicating whether to cluster genes/rows (default: FALSE)
#' @param cluster_cols Logical indicating whether to cluster identities/columns (default: FALSE)
#' @param show_gene_age Logical indicating whether to show gene age as row annotation (default: TRUE)
#' @param show_gene_ids Logical indicating whether to show gene names (default: FALSE)
#' @param gene_annotation Data frame with custom gene annotations, rownames should match gene IDs (default: NULL)
#' @param gene_annotation_colors Named list of color vectors for custom gene annotations (default: NULL)
#' @param annotation_col Data frame with column annotations (default: NULL)
#' @param annotation_col_colors List of colors for column annotations (default: NULL)
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
                                   gene_annotation = NULL,
                                   gene_annotation_colors = NULL,
                                   annotation_col = NULL,
                                   annotation_col_colors = NULL,
                                   ...) {
    
    # Apply log transformation
    e <- expression_matrix |> log1p()
    
    # Filter genes if specific genes are provided
    if (!is.null(genes) && length(genes) > 0) {
        # Check which genes are present in the expression matrix
        present_genes <- intersect(genes, rownames(e))
        missing_genes <- setdiff(genes, rownames(e))
        
        if (length(missing_genes) > 0) {
            warning(sprintf("The following %d gene(s) were not found in the expression matrix: %s", 
                          length(missing_genes), 
                          paste(missing_genes, collapse = ", ")))
        }
        
        if (length(present_genes) == 0) {
            stop("None of the specified genes were found in the expression matrix")
        }
        
        e <- e[present_genes, , drop = FALSE]
        message(sprintf("Using %d out of %d specified genes", length(present_genes), length(genes)))
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
    
    # Handle custom gene annotations first (takes precedence)
    if (!is.null(gene_annotation)) {
        # Filter gene_annotation to only include genes present in the expression matrix
        present_genes <- intersect(rownames(e), rownames(gene_annotation))
        
        if (length(present_genes) > 0) {
            annotation_row <- gene_annotation[present_genes, , drop = FALSE]
            
            # Use provided colors or generate defaults
            if (!is.null(gene_annotation_colors)) {
                annotation_colors <- gene_annotation_colors
            } else {
                # Generate default colors for each annotation column
                annotation_colors <- list()
                for (col_name in colnames(annotation_row)) {
                    unique_vals <- unique(annotation_row[[col_name]])
                    if (is.numeric(unique_vals)) {
                        # For numeric annotations, use a color gradient
                        annotation_colors[[col_name]] <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
                    } else {
                        # For categorical annotations, use distinct colors
                        n_colors <- length(unique_vals)
                        colors <- grDevices::colorRampPalette(c("red","blue","green","orange"))(n_colors)
                        names(colors) <- unique_vals
                        annotation_colors[[col_name]] <- colors
                    }
                }
            }
        } else {
            warning("No genes from gene_annotation found in the expression matrix")
        }
    } else if (show_gene_age) {
        # Use gene age annotation if no custom annotation provided
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
        annotation_col = annotation_col,
        annotation_colors = c(annotation_colors, annotation_col_colors),
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
                                                                 gene_annotation = NULL,
                                                                 gene_annotation_colors = NULL,
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
        gene_annotation = gene_annotation,
        gene_annotation_colors = gene_annotation_colors,
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
                                                               gene_annotation = NULL,
                                                               gene_annotation_colors = NULL,
                                                               ...) {
    
    # Select expression data
    cell_type_annotation <- NULL
    if (show_reps) {
        # Use downsample_expression to get downsampled matrix
        expression_matrix <- downsample_expression(phyex_set@expression, phyex_set@groups, max_cells_per_type)
        
        # Extract cell type information from column names using selected_idents
        cell_names <- colnames(expression_matrix)
        cell_indices <- match(cell_names, colnames(phyex_set@expression))
        cell_types <- as.character(phyex_set@groups[cell_indices])
        names(cell_types) <- cell_names
        
        # Use the current identities label as the column name
        annotation_col_name <- phyex_set@identities_label
        cell_type_annotation <- data.frame(
            row.names = colnames(expression_matrix)
        )
        cell_type_annotation[[annotation_col_name]] <- cell_types[colnames(expression_matrix)]
        
        # Sort expression matrix columns by cell type to group cells together
        cell_type_order <- order(cell_type_annotation[[annotation_col_name]])
        expression_matrix <- expression_matrix[, cell_type_order, drop = FALSE]
        cell_type_annotation <- cell_type_annotation[cell_type_order, , drop = FALSE]
        
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

    # Create column annotation colors if we have cell type annotations
    annotation_col_colors <- NULL
    if (!is.null(cell_type_annotation)) {
        # Get the annotation column name (should match identities_label)
        annotation_col_name <- names(cell_type_annotation)[1]
        
        # Get or create colors for cell types
        cell_type_levels <- unique(cell_type_annotation[[annotation_col_name]])
        if (!is.null(phyex_set@idents_colours) && phyex_set@identities_label %in% names(phyex_set@idents_colours)) {
            # Use existing colors if available
            cell_type_colors <- phyex_set@idents_colours[[phyex_set@identities_label]]
        } else {
            # Generate default colors
            cell_type_colors <- grDevices::colorRampPalette(c("red","blue","green","orange"))(length(cell_type_levels))
            names(cell_type_colors) <- cell_type_levels
        }
        
        # Create the annotation colors list with the dynamic column name
        annotation_col_colors <- list()
        annotation_col_colors[[annotation_col_name]] <- cell_type_colors
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
        gene_annotation = gene_annotation,
        gene_annotation_colors = gene_annotation_colors,
        annotation_col = cell_type_annotation,
        annotation_col_colors = annotation_col_colors,
        ...
    )
}