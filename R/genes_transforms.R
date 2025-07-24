#' @title Gene Expression Transformation Functions
#' @description Collection of functions for transforming gene expression data
#' in PhyloExpressionSet objects.
#' 

#' @title Set Expression Matrix in PhyloExpressionSet
#' @description Generic function to set the expression matrix in a PhyloExpressionSet object.
#' @param phyex_set A PhyloExpressionSet object
#' @param new_expression Matrix to set as the new expression
#' @param ... Additional arguments
#' @return A PhyloExpressionSet object with updated expression
#' @export
set_expression <- S7::new_generic("set_expression", "phyex_set")

#' @export
S7::method(set_expression, BulkPhyloExpressionSet) <- function(phyex_set, new_expression, new_name = NULL, ...) {
    stopifnot(all(dim(new_expression) == dim(phyex_set@.expression)))
    phyex_set@.expression <- new_expression
    if (!is.null(new_name)) phyex_set@name <- new_name
    return(phyex_set)
}

#' @export
S7::method(set_expression, ScPhyloExpressionSet) <- function(phyex_set, new_expression, new_name = NULL, ...) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("Package 'SeuratObject' must be installed to use this function.")
    }
    stopifnot(all(dim(new_expression) == dim(phyex_set@expression)))
    new_seurat <- phyex_set@seurat
    SeuratObject::LayerData(new_seurat, layer = phyex_set@layer) <- new_expression
    phylomap <- data.frame(
        Stratum = as.numeric(phyex_set@strata),
        GeneID = phyex_set@gene_ids
    )
    res <- as_ScPhyloExpressionSet(
        seurat = new_seurat,
        phylomap = phylomap,
        layer = phyex_set@layer,
        name = if (!is.null(new_name)) new_name else phyex_set@name,
        min_cells_per_identity = phyex_set@min_cells_per_identity,
        species = phyex_set@species,
        index_type = phyex_set@index_type
    )
    return(res)
}

#' @title Transform Expression Counts in PhyloExpressionSet
#' @description Apply a transformation function to the expression counts in a PhyloExpressionSet.
#' @param phyex_set A PhyloExpressionSet object
#' @param FUN Function to apply
#' @param FUN_name Name of the transformation function (optional)
#' @param new_name Name for the new dataset (optional)
#' @param ... Additional arguments passed to FUN
#' @return A PhyloExpressionSet object with transformed expression data
#' @export
transform_counts <- function(phyex_set, FUN, FUN_name = deparse(substitute(FUN)), new_name = NULL, ...) {
    f <- match.fun(FUN)
    new_expr <- f(phyex_set@expression, ...)
    if (is.null(new_name)) new_name <- paste(phyex_set@name, "transformed by", FUN_name)
    set_expression(phyex_set, new_expr, new_name = new_name)
}

#' @title Short Alias for Transform Counts
#' @description Convenience alias for transform_counts function.
#' @export
tf <- transform_counts

#' @title Normalise Stage Expression Data
#' @description Normalise expression data to a specified total expression level per sample.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param total Numeric value to normalise each sample to (default: 1e6)
#' 
#' @return A PhyloExpressionSet object with normalised expression data
#' 
#' @examples
#' # Normalise to 1 million total expression per sample
#' # normalised_set <- normalise_stage_expression(phyex_set, total = 1e6)
#' 
#' @export
normalise_stage_expression <- function(phyex_set, total=1e6) {
    new_expr <- sweep(phyex_set@expression, 2, colSums(phyex_set@expression), FUN="/") * total
    set_expression(phyex_set, new_expr, new_name = paste(phyex_set@name, "normalised to", total))
}

#' @title Standardise Expression Data
#' @description Standardise gene expression data by centring and scaling each gene.
#' 
#' @param e Matrix of expression values with genes as rows and samples as columns
#' 
#' @return Matrix with standardised expression values (mean=0, sd=1 for each gene)
#' 
#' @details
#' This function standardises each gene's expression profile by subtracting the mean
#' and dividing by the standard deviation. Genes with zero or undefined variance
#' are set to zero. This is useful for comparing expression patterns across genes
#' with different absolute expression levels.
#' 
#' @examples
#' # Standardise expression data
#' # std_expr <- .to_std_expr(expression_matrix)
#' 
#' @keywords internal
.to_std_expr <- function(e) {
    row_sd <- apply(e, 1, stats::sd, na.rm = TRUE)
    valid <- row_sd > 0 & is.finite(row_sd)
    e[valid, ] <- t(scale(t(e[valid, , drop = FALSE]), center = TRUE, scale = TRUE))
    e[!valid, ] <- 0
    e
}
