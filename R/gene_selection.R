
# expressed genes

#' @export
top_variance_genes <- function(phyex_set, p = .99){
    avg_counts <- rowVars(phyex_set@counts)
    names(avg_counts) <- phyex_set@gene_ids
    
    top_genes <- names(avg_counts)[avg_counts >= stats::quantile(avg_counts, p)]
    
    return(top_genes)
}

#' @export
top_expression_genes <- function(phyex_set, p = .99){
    avg_counts <- rowMeans(phyex_set@counts)
    names(avg_counts) <- phyex_set@gene_ids
    
    
    top_genes <- names(avg_counts)[avg_counts >= stats::quantile(avg_counts, p)]

    return(top_genes)
}