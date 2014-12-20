#' @title Compute a permutation matrix for test statistics
#' @description This function computes the TAI for a row permutated PhyloExpressionSet or DivergenceExpressionSet.
#'
#' One can specify the number of permutations which corresponds to the number of \code{\link{TAI}} or \code{\link{TDI}} profiles
#' that are being returned as data matrix. The function then returns a \code{TAI} or \code{TDI} matrix holding
#' the \code{TAI} or \code{TDI} profiles of the permutated PhyloExpressionSets or DivergenceExpressionSets. This procedure
#' can be used for building test statistics based on the \code{TAI} or \code{TDI} profiles. 
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param permutations a numeric value specifying the number of permutations to be performed.
#' @details The sampled \code{\link{TAI}} or \code{\link{TDI}} matrix samples the phylostratum or divergence-stratum vector of
#' a given PhyloExpressionSet or DivergenceExpressionSet and computes the corresponding \code{TAI} or \code{TDI} profiles
#' of the randomly assigned phylostrata or divergence-strata. This sampling is then performed N times, yielding N randomly sampled \code{TAI} or \code{TDI} profiles.
#' This random \code{TAI} or \code{TDI} profile matrix can then be used to perform statistical tests 
#' (such as the \code{\link{FlatLineTest}} or \code{\link{ReductiveHourglassTest}}) based on the significance of \code{TAI} or \code{TDI} patterns.
#' @return a numeric matrix representing N randomly permuted \code{TAI} or \code{TDI} profiles.
#' @references Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}}
#' @examples 
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#' 
#' # example PhyloExpressionSet using 100 permutations
#' randomTAI.Matrix <- bootMatrix(PhyloExpressionSetExample, permutations = 100)
#' 
#' # example DivergenceExpressionSet using 100 permutations
#' randomTDI.Matrix <- bootMatrix(DivergenceExpressionSetExample, permutations = 100)
#' 
#' 
#' @export
bootMatrix <- function(ExpressionSet,permutations = 1000)
{
        
        is.ExpressionSet(ExpressionSet)
        
        nCols <- dim(ExpressionSet)[2]
        bootstrapMatrix <- matrix(NA_real_, permutations, (nCols - 2))
        ExprSet <- as.matrix(ExpressionSet[ , 3:nCols])
        AgeVector <- as.vector(ExpressionSet[ , 1])
        
        bootstrapMatrix <- cpp_bootMatrix(ExprSet, AgeVector, permutations)
        colnames(bootstrapMatrix) <- names(ExpressionSet)[3:nCols]
        rownames(bootstrapMatrix) <- 1:permutations
        
        return(bootstrapMatrix)
        
        
}