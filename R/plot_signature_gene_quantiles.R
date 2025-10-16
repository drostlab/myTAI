
#' @title Plot Signature Across Gene Expression Quantiles
#' @description Create a plot showing how the transcriptomic signature changes when
#' genes are progressively removed based on expression quantiles.
#' 
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)
#' @param quantiles Numeric vector of quantiles to test (default: c(1.0, 0.99, 0.95, 0.90, 0.80))
#' @param selection_FUN Function to select genes for removal (default: genes_top_mean)
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
#' The analysis works with both bulk and single-cell data, helping to determine
#' whether phylotranscriptomic patterns are driven by a few highly expressed genes
#' or represent broad transcriptomic trends.
#' 
#' @examples
#' # Plot signature across expression quantiles for bulk data
#' phyex_set <- example_phyex_set |>
#'     select_genes(example_phyex_set@gene_ids[1:100])
#' phyex_set@null_conservation_sample_size <- 500
#' p <- plot_signature_gene_quantiles(phyex_set, quantiles = c(0.95, 0.90))
#' 
#' @import purrr
#' @export
plot_signature_gene_quantiles <- function(phyex_set,
                                          quantiles=c(1.0, 0.99, 0.95, 0.90, 0.80),
                                          selection_FUN=genes_top_mean,
                                          ...) {
    phyex_sets <- quantiles |>
        map(\(p) selection_FUN(phyex_set, p=p)) |>
        map(\(genes) remove_genes(phyex_set, genes=genes)) |>
        map2(quantiles, \(set, p) {set@name <- as.character(p); set})
    
    plot_signature_multiple(phyex_sets, legend_title="Quantile", ...)
}