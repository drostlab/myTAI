
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