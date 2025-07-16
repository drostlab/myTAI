
#' @title Plot Sample Space Visualization
#' @description Create a dimensional reduction plot to visualize sample relationships
#' in gene expression space using PCA or UMAP.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param method Character string specifying the dimensionality reduction method: 
#' "PCA" or "UMAP" (default: "PCA")
#' 
#' @return A ggplot2 object showing the sample space visualization
#' 
#' @details
#' This function performs log1p transformation on expression data, removes genes with
#' zero variance, and applies the specified dimensionality reduction method. Samples
#' are colored by their group assignments.
#' 
#' @examples
#' # Create PCA plot
#' # pca_plot <- plot_sample_space(phyex_set, method = "PCA")
#' 
#' # Create UMAP plot
#' # umap_plot <- plot_sample_space(phyex_set, method = "UMAP")
#' 
#' @import ggplot2
#' @importFrom uwot umap
#' @export
plot_sample_space <- function(phyex_set, method=c("PCA", "UMAP")) {
    method <- match.arg(method)
    
    expr <- log1p(phyex_set@counts)
    # Remove genes with zero variance
    nonzero_var_genes <- apply(expr, 1, function(x) var(x) != 0)
    expr <- expr[nonzero_var_genes, , drop = FALSE]
    
    if (method == "PCA") {
        coords <- prcomp(t(expr), scale. = TRUE)$x[, 1:2]
    }
    else if (method == "UMAP"){
        coords <- uwot::umap(t(expr), scale = TRUE)
    }
    
    df <- as.data.frame(coords)
    colnames(df) <- c("V1", "V2")
    df$Group <- phyex_set@groups
    
    ggplot(df, aes(V1, V2, colour=Group)) +
        geom_point(size=3) +
        ggtitle(paste("Sample", method)) + 
        scale_colour_viridis_d() +
        theme_minimal()
}