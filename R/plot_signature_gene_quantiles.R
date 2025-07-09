
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