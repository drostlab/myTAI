#' @title Gene Expression Filtering Functions
#' @description Collection of functions for filtering genes based on expression patterns
#' in PhyloExpressionSet objects.
#' 

#' @title Select Top Genes by Expression Metric
#' @description Generic function to select genes with the highest values for a given expression metric.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param FUN Function to calculate gene-wise expression metric (default: rowMeans)
#' @param p Quantile threshold for gene selection (default: 0.99)
#' @param ... Additional arguments passed to FUN
#' 
#' @return Character vector of gene IDs with metric values >= p quantile
#' 
#' @details
#' This function applies the specified function to calculate a metric for each gene
#' across samples, then selects genes above the specified quantile threshold.
#' Common functions include rowMeans for mean expression, rowVars for variance,
#' or custom functions for other metrics.
#' 
#' @examples
#' # Select top 1% most expressed genes by mean
#' high_expr_genes <- genes_top_expr(example_phyex_set, function(x) apply(x, 1, mean), p = 0.99)
#' 
#' # Select top 5% most variable genes
#' high_var_genes <- genes_top_expr(example_phyex_set, function(x) apply(x, 1, var), p = 0.95)
#' 
#' # Select top genes by median expression
#' high_median_genes <- genes_top_expr(example_phyex_set, function(x) apply(x, 1, median), p = 0.9)
#' 
#' @export
genes_top_expr <- function(phyex_set, FUN = rowMeans, p = .99, ...) {
    check_PhyloExpressionSet(phyex_set)
    FUN <- match.fun(FUN)
    metric_values <- FUN(phyex_set@expression_collapsed, ...)
    names(metric_values) <- phyex_set@gene_ids

    if (p >= 1) {
        return(character(0)) # Return no genes if p = 1
    }
    if (p <= 0) {
        return(names(metric_values)) # Return all genes if p = 0
    }
    
    top_genes <- names(metric_values)[metric_values > stats::quantile(metric_values, p, na.rm = TRUE)]
    
    return(top_genes)
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
#' high_var_genes <- genes_top_variance(example_phyex_set, p = 0.99)
#' 
#' @export
genes_top_variance <- function(phyex_set, p = .99) {
    genes_top_expr(phyex_set, FUN = rowVars, p = p)
}

#' @title Select Top Mean Expressed Genes
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
#' # Select top 1% most expressed genes by mean
#' high_expr_genes <- genes_top_mean(example_phyex_set, p = 0.99)
#' 
#' @export
genes_top_mean <- function(phyex_set, p = .99) {
    genes_top_expr(phyex_set, FUN = rowMeans, p = p)
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
#' low_expr_genes <- genes_lowly_expressed(example_phyex_set, threshold = 1)
#' 
#' @export
genes_lowly_expressed <- function(phyex_set, threshold = 1) {
    avg_counts <- rowMeans(phyex_set@expression)
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
#' # Filter top 10% most variable genes
#' # filtered_expr <- genes_filter_dynamic(expression_matrix, thr = 0.9)
#' 
#' @keywords internal
genes_filter_dynamic <- function(e, thr=0.9) {
    var_genes <- apply(e, 1, stats::var, na.rm = TRUE)
    var_genes <- var_genes[!is.na(var_genes) & !is.nan(var_genes)]
    if (length(var_genes) == 0) {
        return(e[FALSE, , drop = FALSE])  # Return empty matrix with same structure
    }
    cutoff <- stats::quantile(var_genes, thr, na.rm = TRUE)
    all_var_genes <- apply(e, 1, stats::var, na.rm = TRUE)
    valid_genes <- !is.na(all_var_genes) & !is.nan(all_var_genes) & all_var_genes > cutoff
    e[valid_genes, , drop = FALSE]
}
