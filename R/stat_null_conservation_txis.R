
#' @title Generate Null Conservation TXI Distribution
#' @description Generate a null distribution of transcriptomic index values for
#' conservation testing by permuting phylostratum assignments.
#' 
#' @param strata_vector Numeric vector of phylostratum assignments
#' @param count_matrix Matrix of expression counts
#' @param sample_size Number of permutations to generate
#' 
#' @return Matrix of permuted TXI values for null hypothesis testing
#' 
#' @details
#' This function creates a null distribution by randomly permuting phylostratum
#' assignments while keeping the expression matrix fixed. This preserves the
#' expression structure while breaking the phylostratum-expression relationship.
#' # Generate null TXI distribution
#' # null_txis <- stat_generate_conservation_txis(strata_vec, expr_matrix, 1000)
#' 
#' @keywords internal
stat_generate_conservation_txis <- function(strata_vector,
                                            count_matrix,
                                            sample_size) {
    
    permuted_txi_matrix <- cpp_nullTXIs(count_matrix,
                                        strata_vector,
                                        sample_size)
    
    colnames(permuted_txi_matrix) <- colnames(count_matrix)
    return(permuted_txi_matrix)
}

#' @title Memoized Null Conservation TXI Generation
#' @description Memoized version of generate_conservation_txis for improved performance
#' with repeated calls using the same parameters.
#' 
#' @details
#' This function caches results to avoid recomputing expensive permutations
#' when the same parameters are used multiple times.
#' 
#' @keywords internal
#' @importFrom memoise memoise
.memo_generate_conservation_txis <- memoise::memoise(stat_generate_conservation_txis)