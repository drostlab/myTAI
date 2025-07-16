
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

#' @title Standardize Expression Data
#' @description Standardize gene expression data by centering and scaling each gene.
#' 
#' @param e Matrix of expression values with genes as rows and samples as columns
#' 
#' @return Matrix with standardized expression values (mean=0, sd=1 for each gene)
#' 
#' @details
#' This function standardizes each gene's expression profile by subtracting the mean
#' and dividing by the standard deviation. Genes with zero or undefined variance
#' are set to zero. This is useful for comparing expression patterns across genes
#' with different absolute expression levels.
#' 
#' @examples
#' # Standardize expression data
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