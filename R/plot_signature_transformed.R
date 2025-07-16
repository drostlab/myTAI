
#' @title Count Transformation Functions
#' @description Predefined list of transformation functions for count data normalization.
#' 
#' @format A named list of transformation functions:
#' \describe{
#'   \item{none}{Identity transformation (no change)}
#'   \item{sqrt}{Square root transformation}
#'   \item{log2}{log2(x+1) transformation}
#'   \item{vst}{Variance stabilizing transformation (DESeq2)}
#'   \item{rlog}{Regularized log transformation (DESeq2)}
#'   \item{rank}{Rank transformation within each sample}
#' }
#' @export
COUNT_TRANSFORMS <- list(
    "none" = identity,
    "sqrt" = sqrt,
    "log2" = \(x) log2(x+1),
    "vst" = \(x) DESeq2::vst(round(x, digits=0)),
    "rlog" = \(x) DESeq2::rlogTransformation(round(x, digits=0)),
    "rank" = \(x) apply(x, 2, base::rank)
)

#' @title Plot Signature Under Different Transformations
#' @description Compare transcriptomic signatures under various data transformations
#' to assess the robustness of phylotranscriptomic patterns.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param transformations Named list of transformation functions (default: COUNT_TRANSFORMS)
#' @param ... Additional arguments passed to plot_signature_multiple
#' 
#' @return A ggplot2 object showing signatures under different transformations
#' 
#' @details
#' This function applies different transformations to the same dataset and
#' compares the resulting transcriptomic signatures. This is useful for assessing
#' whether phylotranscriptomic patterns are robust to different data processing
#' approaches or are artifacts of specific transformations.
#' 
#' @examples
#' # Compare signatures under all default transformations
#' # p <- plot_signature_transformed(phyex_set)
#' 
#' # Custom transformations
#' # custom_transforms <- list(raw = identity, log = log1p)
#' # p2 <- plot_signature_transformed(phyex_set, transformations = custom_transforms)
#' 
#' @export
plot_signature_transformed <- function(phyex_set,
                                       transformations=COUNT_TRANSFORMS,
                                       ...) {
    phyex_sets <- transformations |>
        purrr::imap(\(x, idx) transform_counts(phyex_set, FUN=x, FUN_name=idx))
    
    plot_signature_multiple(phyex_sets, legend_title="Transformation", ...)
}