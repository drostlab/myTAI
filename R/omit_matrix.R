#' @title Compute TXI Profiles Omitting Each Gene
#' @description For each gene i, exclude the corresponding gene i from the
#' PhyloExpressionSet and compute the TXI profile for the dataset with gene i excluded.
#'  
#' This procedure results in a TXI profile matrix storing the TXI profile for each omitted gene i.
#' 
#' @param phyex_set A PhyloExpressionSet object        
#' @return A numeric matrix storing TXI profiles for each omitted gene i
#' 
#' @details
#' This function systematically removes each gene and recalculates the transcriptomic
#' index profile to assess the contribution of individual genes to the overall pattern.
#' This is useful for identifying genes that have a large influence on the
#' phylotranscriptomic signature.
#' 
#' @examples
#' # Compute omit matrix for a PhyloExpressionSet
#' # omit_mat <- omit_matrix(phyex_set)
#' 
#' @author Hajk-Georg Drost
#' @export
omit_matrix <- function(phyex_set) {

    oMatrix <- .omit_matrix(phyex_set@counts_collapsed, phyex_set@strata)
    
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
