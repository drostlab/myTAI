#' @title Gene Expression Transformation Functions
#' @description Collection of functions for transforming and filtering gene expression data
#' in PhyloExpressionSet objects.
#' 
#' @importFrom matrixStats rowVars

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
S7::method(transform_counts, PhyloExpressionSet) <- function(phyex_set, 
                                                             FUN,
                                                             FUN_name=deparse(substitute(FUN)),
                                                             new_name=paste(phyex_set@name, "transformed by", FUN_name),
                                                             ...) {
    f <- match.fun(FUN)
    phyex_set@counts <- f(phyex_set@counts, ...)
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
    phyex_set@counts <- sweep(phyex_set@counts, 2, colSums(phyex_set@counts), FUN="/") * total
    phyex_set
}

#' @title Select Top Variable Genes
#' @description Select genes with the highest variance across samples.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param p Quantile threshold for gene selection (default: 0.99)
#' 
#' @return Character vector of gene IDs with variance >= p quantile
#' 
#' @details
#' This function identifies genes with the highest variance across samples,
#' which are often the most informative for downstream analyses.
#' 
#' @examples
#' # Select top 1% most variable genes
#' # high_var_genes <- top_variance_genes(phyex_set, p = 0.99)
#' 
#' @export
top_variance_genes <- function(phyex_set, p = .99) {
    avg_counts <- rowVars(phyex_set@counts)
    names(avg_counts) <- phyex_set@gene_ids
    
    top_genes <- names(avg_counts)[avg_counts >= stats::quantile(avg_counts, p)]
    
    return(top_genes)
}

#' @title Select Top Expressed Genes
#' @description Select genes with the highest mean expression across samples.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param p Quantile threshold for gene selection (default: 0.99)
#' 
#' @return Character vector of gene IDs with mean expression >= p quantile
#' 
#' @details
#' This function identifies genes with the highest mean expression levels,
#' which are often the most reliably detected and functionally important.
#' 
#' @examples
#' # Select top 1% most expressed genes
#' # high_expr_genes <- top_expression_genes(phyex_set, p = 0.99)
#' 
#' @export
top_expression_genes <- function(phyex_set, p = .99) {
    avg_counts <- rowMeans(phyex_set@counts)
    names(avg_counts) <- phyex_set@gene_ids
    
    top_genes <- names(avg_counts)[avg_counts >= stats::quantile(avg_counts, p)]

    return(top_genes)
}

#' @title Select Lowly Expressed Genes
#' @description Select genes with mean expression below a specified threshold.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param threshold Mean expression threshold (default: 1)
#' 
#' @return Character vector of gene IDs with mean expression <= threshold
#' 
#' @details
#' This function identifies genes with low mean expression levels, which
#' might be candidates for filtering or separate analysis.
#' 
#' @examples
#' # Select genes with mean expression <= 1
#' # low_expr_genes <- lowly_expressed_genes(phyex_set, threshold = 1)
#' 
#' @export
lowly_expressed_genes <- function(phyex_set, threshold = 1) {
    avg_counts <- rowMeans(phyex_set@counts)
    names(avg_counts) <- phyex_set@gene_ids
    
    genes <- names(avg_counts)[avg_counts <= threshold]
    
    return(genes)
}

#' @title Filter Dynamic Expression Genes
#' @description Filter genes based on expression variance to select the most dynamically expressed genes.
#' 
#' @param e Matrix of expression values with genes as rows and samples as columns
#' @param thr Threshold quantile for variance filtering (default: 0.9)
#' 
#' @return Matrix containing only genes above the variance threshold
#' 
#' @details
#' This function calculates the variance for each gene across samples and retains
#' only genes with variance above the specified quantile threshold. This helps
#' focus analysis on genes that show significant expression changes.
#' 
#' @examples
#' # Filter top 10% most variable genes
#' # filtered_expr <- filter_dyn_expr(expression_matrix, thr = 0.9)
#' 
#' @keywords internal
filter_dyn_expr <- function(e, thr=0.9) {
    var_genes <- apply(e, 1, stats::var)
    cutoff <- stats::quantile(var_genes, thr)
    e[var_genes > cutoff, ]
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
