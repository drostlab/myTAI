#' @title Compute TAI or TDI Profiles Omitting a Given Gene
#' @description For each gene i, exclude the corresponding gene i from the global
#'  PhyloExpressionSet or DivergenceExpressionSet and compute the \code{\link{TAI}} or \code{\link{TDI}} 
#'  profile for the corresponding global PhyloExpressionSet or DivergenceExpressionSet
#'  with excluded gene i. 
#'  
#'  This procedure results in a TAI or TDI profile Matrix storing the TAI or TDI profile for each omitted gene i.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.        
#' @return a numeric matrix storing TAI or TDI profile for each omitted gene i.
#' @author Hajk-Georg Drost
#' @examplesIf FALSE
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#'
#' # example PhyloExpressionSet
#' omMatrix_ps <- omitMatrix(PhyloExpressionSetExample)
#'
#' # example DivergenceExpressionSet
#' omMatrix_ds <- omitMatrix(DivergenceExpressionSetExample)
#' @export
omit_matrix <- function(phyex_set) {

    oMatrix <- .omit_matrix(phyex_set@count_matrix, phyex_set@strata_vector)
    
    colnames(oMatrix) <- phyex_set@conditions
    rownames(oMatrix) <- paste0("(-) ", phyex_set@gene_ids)
    
    return(oMatrix)
        
}

.omit_matrix <- function(count_matrix, strata_vector) {
    n <- nrow(count_matrix)
    m <- ncol(count_matrix)
    numerator <- matrix(t(strata_vector) %*% count_matrix, n, m, byrow = TRUE) - count_matrix * strata_vector
    denominator <- matrix(colSums(count_matrix), n, m, byrow = TRUE) - count_matrix
    
    return(numerator / denominator)
}
