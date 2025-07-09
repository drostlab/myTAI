
COUNT_TRANSFORMS <- list(
    "none" = identity,
    "sqrt" = sqrt,
    "log2" = \(x) log2(x+1),
    "vst" = \(x) DESeq2::vst(round(x, digits=0)),
    "rlog" = \ (x) DESeq2::rlogTransformation(round(x, digits=0)),
    "rank" = \(x) apply(x, 2, base::rank)
)

#' @export
plot_signature_transformed <- function(phyex_set,
                                       transformations=COUNT_TRANSFORMS,
                                       ...) {
    phyex_sets <- transformations |>
        purrr::imap(\(x, idx) transform_counts(phyex_set, FUN=x, FUN_name=idx))
    
    plot_signature_multiple(phyex_sets, legend_title="Transformation", ...)
}