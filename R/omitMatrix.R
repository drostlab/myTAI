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
#' @examples
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
#' 
#' 
#' @export
omitMatrix <- function(ExpressionSet)
{
        
        is.ExpressionSet(ExpressionSet)
        
        ncols <- dim(ExpressionSet)[2]
        
        oMatrix <- matrix(NA_real_, ncol = (ncols - 2), nrow = nrow(ExpressionSet))
        oMatrix <- cpp_omitMatrix(as.matrix(ExpressionSet[ , 3:ncols]),as.vector(ExpressionSet[ , 1]))
        
        colnames(oMatrix) <- names(ExpressionSet)[3:ncols]
        rownames(oMatrix) <- paste0("(-) ",ExpressionSet[ , 2])
        
        return(oMatrix)
        
}