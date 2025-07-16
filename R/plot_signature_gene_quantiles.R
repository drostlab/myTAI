
#' @title Plot Signature Across Gene Expression Quantiles
#' @description Create a plot showing how the transcriptomic signature changes when
#' genes are progressively removed based on expression quantiles.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param quantiles Numeric vector of quantiles to test (default: c(1.0, 0.99, 0.95, 0.90, 0.80))
#' @param selection_FUN Function to select genes for removal (default: top_expression_genes)
#' @param ... Additional arguments passed to plot_signature_multiple
#' 
#' @return A ggplot2 object showing signatures across different quantiles
#' 
#' @details
#' This function systematically removes genes based on expression quantiles and
#' shows how the transcriptomic signature changes. This is useful for understanding
#' the contribution of highly expressed genes to the overall pattern and for
#' assessing the robustness of phylotranscriptomic patterns.
#' 
#' @examples
#' # Plot signature across expression quantiles
#' # p <- plot_signature_gene_quantiles(phyex_set)
#' 
#' # Custom quantiles
#' # p2 <- plot_signature_gene_quantiles(phyex_set, quantiles = c(1.0, 0.95, 0.90))
#' 
#' @import purrr
#' @export
plot_signature_gene_quantiles <- function(phyex_set,
                                          quantiles=c(1.0, 0.99, 0.95, 0.90, 0.80),
                                          selection_FUN=top_expression_genes,
                                          ...) {
    phyex_sets <- quantiles |>
        map(\(p) selection_FUN(phyex_set, p=p)) |>
        map(\(genes) remove_genes(phyex_set, genes=genes)) |>
        map2(quantiles, \(set, p) {set@name <- as.character(p); set})
    
    plot_signature_multiple(phyex_sets, legend_title="Quantile", ...)
}