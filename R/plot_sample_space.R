
#' @title Plot Sample Space Visualization
#' @description Create a dimensional reduction plot to visualize sample relationships
#' in gene expression space using PCA or UMAP.
#' 
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)
#' @param method Character string specifying the dimensionality reduction method: 
#' "PCA" or "UMAP" (default: "PCA")
#' @param colour_by Character string specifying what to colour by: "identity" (default), 
#' "TXI"
#' @param seed Integer seed for reproducible UMAP results (default: 42)
#' @param ... Additional arguments passed to specific methods
#' 
#' @return A ggplot2 object showing the sample space visualisation
#' 
#' @details
#' This function performs log1p transformation on expression data, removes genes with
#' zero variance, and applies the specified dimensionality reduction method. Samples
#' are coloured by their group assignments or TAI values.
#' 
#' @examples
#' # Create PCA plot coloured by identity
#' pca_plot <- plot_sample_space(example_phyex_set, method = "PCA", colour_by = "identity")
#' 
#' # Create UMAP plot coloured by TXI
#' umap_plot <- plot_sample_space(example_phyex_set, method = "UMAP", colour_by = "TXI")
#' 
#' @import ggplot2
#' @export
plot_sample_space <- S7::new_generic("plot_sample_space", "phyex_set",
    function(phyex_set, 
             method = c("PCA", "UMAP"), 
             colour_by = c("identity", "TXI"), 
             seed = 42,
             ...) {
        S7::S7_dispatch()
    }
)

#' @export
S7::method(plot_sample_space, BulkPhyloExpressionSet) <- function(phyex_set, 
                                                                method = c("PCA", "UMAP"), 
                                                                colour_by = c("identity", "TXI"), 
                                                                seed = 42,
                                                                ...) {
    method <- match.arg(method)
    colour_by <- match.arg(colour_by)
    
    expr <- log1p(phyex_set@expression)
    # Remove genes with zero variance
    nonzero_var_genes <- apply(expr, 1, function(x) stats::var(x) != 0)
    expr <- expr[nonzero_var_genes, , drop = FALSE]
    
    if (method == "PCA") {
        coords <- stats::prcomp(t(expr), scale. = TRUE)$x[, 1:2]
    }
    else if (method == "UMAP"){
        if (!requireNamespace("uwot", quietly = TRUE)) {
            stop("Package 'uwot' must be installed to compute UMAP.")
        }
        set.seed(seed)
        coords <- uwot::umap(t(expr), scale = TRUE)
    }
    
    df <- as.data.frame(coords)
    colnames(df) <- c("V1", "V2")
    
    # Set up colouring
    if (colour_by == "identity") {
        df$Colour_Variable <- phyex_set@groups
        colour_label <- phyex_set@identities_label
        use_discrete <- TRUE
    } else if (colour_by == "TXI") {
        df$Colour_Variable <- phyex_set@TXI_sample
        colour_label <- phyex_set@index_full_name
        use_discrete <- FALSE
    }
    
    p <- ggplot(df, aes(V1, V2, colour = Colour_Variable)) +
        geom_point(size = 3) +
        labs(
            title = paste(method, "Plot -", phyex_set@name),
            x = paste(method, "1"),
            y = paste(method, "2"),
            colour = colour_label
        ) +
        theme_minimal()
    
    # Apply appropriate colour scale
    if (use_discrete) {
        p <- p + scale_colour_viridis_d()
    } else {
        p <- p + scale_colour_gradient2(low = "blue", mid = "white", high = "red", 
                                       midpoint = median(df$Colour_Variable, na.rm = TRUE))
    }
    
    return(p)
}

#' @export
S7::method(plot_sample_space, ScPhyloExpressionSet) <- function(phyex_set, 
                                                               method = c("PCA", "UMAP"), 
                                                               colour_by = c("identity", "TXI"), 
                                                               seed = 42,
                                                               ...) {
    method <- match.arg(method)
    colour_by <- match.arg(colour_by)
    
    # Check if reduction is available in stored reductions
    reduction_key <- tolower(method)
    coords <- NULL
    
    if (reduction_key %in% names(phyex_set@reductions)) {
        # Use existing reduction
        coords <- phyex_set@reductions[[reduction_key]][, 1:2, drop = FALSE]
    } else {
        # Compute reduction if not available
        coords <- .compute_reduction(phyex_set@expression, method = method, seed = seed)
    }
    
    # Prepare plotting data
    plot_data <- data.frame(
        V1 = coords[, 1],
        V2 = coords[, 2],
        cell_id = rownames(coords)
    )
    
    # Set up colouring
    if (colour_by == "identity") {
        plot_data$Colour_Variable <- phyex_set@groups[match(plot_data$cell_id, names(phyex_set@groups))]
        colour_label <- phyex_set@identities_label
        use_discrete <- TRUE
    } else if (colour_by == "TXI") {
        plot_data$Colour_Variable <- phyex_set@TXI_sample[match(plot_data$cell_id, names(phyex_set@TXI_sample))]
        colour_label <- phyex_set@index_full_name
        use_discrete <- FALSE
    }
    
    # Filter out cells with NA values
    plot_data <- plot_data[!is.na(plot_data$Colour_Variable), ]
    
    p <- ggplot(plot_data, aes(V1, V2, colour = Colour_Variable)) +
        geom_point(size = 0.5, alpha = 0.8) +
        labs(
            title = paste(method, "Plot -", phyex_set@name),
            x = paste(method, "1"),
            y = paste(method, "2"),
            colour = colour_label
        ) +
        theme_minimal() +
        theme(
            legend.position = "right",
            panel.grid = element_blank()
        )
        
    # Apply appropriate colour scale
    if (use_discrete) {
        p <- p + scale_colour_viridis_d()
    } else {
        p <- p + scale_colour_gradient2(low = "blue", mid = "white", high = "red", 
                                        midpoint = median(plot_data$Colour_Variable, na.rm = TRUE))
    }
    
    return(p)
}