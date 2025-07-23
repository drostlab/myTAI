#' @title Gene Expression Transformation Functions
#' @description Collection of functions for transforming gene expression data
#' in PhyloExpressionSet objects.
#' 
#' @title Transform Expression Counts in PhyloExpressionSet
#' @description Apply a transformation function to the expression counts in a PhyloExpressionSet.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param ... Additional arguments. For the PhyloExpressionSet method, arguments are FUN (function to apply), FUN_name (character string naming the transformation function), new_name (character string for the new dataset name), and additional arguments passed to the transformation function
#' 
#' @return A PhyloExpressionSet object with transformed expression data
#' 
#' @examples
#' # Apply log transformation
#' # log_set <- transform_counts(phyex_set, log1p, "log1p")
#' 
#' @export
transform_counts <- S7::new_generic("transform_counts", "phyex_set")

#' @export
S7::method(transform_counts, BulkPhyloExpressionSet) <- function(phyex_set, 
                                                                 FUN,
                                                                 FUN_name=deparse(substitute(FUN)),
                                                                 new_name=paste(phyex_set@name, "transformed by", FUN_name),
                                                                 ...) {
    f <- match.fun(FUN)
    phyex_set@expression <- f(phyex_set@expression, ...)
    phyex_set@name <- new_name
    return(phyex_set)
}

#' @export
S7::method(transform_counts, ScPhyloExpressionSet) <- function(phyex_set, 
                                                              FUN, 
                                                              FUN_name = deparse(substitute(FUN)),
                                                              new_name = paste(phyex_set@name, FUN_name),
                                                              ...) {
    
    # Transform the Seurat object
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' must be installed to use this function.")
    }
    
    # Get current expression data
    current_data <- tryCatch({
        if (phyex_set@slot == "data") {
            Seurat::GetAssayData(phyex_set@seurat, layer = "data")
        } else {
            Seurat::GetAssayData(phyex_set@seurat, layer = phyex_set@slot)
        }
    }, error = function(e) {
        Seurat::GetAssayData(phyex_set@seurat, slot = phyex_set@slot)
    })
    
    # Apply transformation
    transformed_data <- FUN(current_data, ...)
    
    # Create new Seurat object with transformed data
    new_seurat <- phyex_set@seurat
    
    # Set the transformed data back
    tryCatch({
        if (phyex_set@slot == "data") {
            new_seurat <- Seurat::SetAssayData(new_seurat, layer = "data", new.data = transformed_data)
        } else {
            new_seurat <- Seurat::SetAssayData(new_seurat, layer = phyex_set@slot, new.data = transformed_data)
        }
    }, error = function(e) {
        new_seurat <- Seurat::SetAssayData(new_seurat, slot = phyex_set@slot, new.data = transformed_data)
    })
    
    # Create phylomap from original object
    phylomap <- data.frame(
        Stratum = as.numeric(phyex_set@strata),
        GeneID = phyex_set@gene_ids
    )
    
    # Create new ScPhyloExpressionSet
    return(as_ScPhyloExpressionSet(
        seurat = new_seurat,
        phylomap = phylomap,
        cell_identity = phyex_set@cell_identity,
        slot = phyex_set@slot,
        name = new_name,
        min_cells_per_type = phyex_set@min_cells_per_type,
        species = phyex_set@species,
        index_type = phyex_set@index_type
    ))
}

#' @title Short Alias for Transform Counts
#' @description Convenience alias for transform_counts function.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param ... Arguments passed to transform_counts
#' 
#' @return A PhyloExpressionSet object with transformed expression data
#' 
#' @examples
#' # Short alias for transformation
#' # log_set <- tf(phyex_set, log1p)
#' 
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
    phyex_set@expression <- sweep(phyex_set@expression, 2, colSums(phyex_set@expression), FUN="/") * total
    phyex_set
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
#' # std_expr <- to_std_expr(expression_matrix)
#' 
#' @keywords internal
to_std_expr <- function(e) {
    row_sd <- apply(e, 1, stats::sd, na.rm = TRUE)
    valid <- row_sd > 0 & is.finite(row_sd)
    e[valid, ] <- t(scale(t(e[valid, , drop = FALSE]), center = TRUE, scale = TRUE))
    e[!valid, ] <- 0
    e
}
