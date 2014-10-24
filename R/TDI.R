#' @title Function to compute the Transcriptome Divergence Index (TDI)
#' @description This function computes the sequence distance based transcriptome divergence index (TDI) introduced by
#' Quint et al., 2012.
#' @param DivergenceExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @details
#' 
#' The TDI measure represents the weighted arithmetic mean (expression levels as
#' weights for the divergence-stratum value) over all gene divergence categories denoted as \emph{divergence-strata}.
#'
#'
#' \deqn{TDI_s = \sum (e_is * ds_i) / \sum e_is}
#' 
#' where TDI_s denotes the TDI value in developmental stage s, e_is denotes the gene expression level of gene i in stage s, and ds_i denotes the corresponding divergence-stratum of gene i, \eqn{i = 1,...,N} and N = total number of genes.
#' 
#' Internally the function is written in C++ to speed up TDI computations.
#' @return a numeric vector containing the TDI values for all given developmental stages.
#' @references 
#' Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{TAI}},  \code{\link{PlotPattern}}, \code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}}
#' @examples \dontrun{
#' 
#' # reading a standard DivergenceExpressionSet
#' data(DivergenceExpressionSetExample)
#'
#' # computing the TDI profile of a given DivergenceExpressionSet object
#' TDIs <- TDI(DivergenceExpressionSetExample)
#' 
#' # you can then visualize the TDI results using PlotPattern()
#' PlotPattern(DivergenceExpressionSetExample, type = "l",
#'             xlab = "Ontogeny", ylab = "TDI", lwd = 9)
#'
#' }
#' @export

TDI <- function(DivergenceExpressionSet)
{
        
        is.ExpressionSet(DivergenceExpressionSet)
        
        nCols <- dim(DivergenceExpressionSet)[2]
        ExpressionMatrix <- DivergenceExpressionSet[ , 3:nCols]
        Divergencestratum <- DivergenceExpressionSet[ , 1]
        TDIProfile <- vector(mode = "numeric",length = nCols-2)
        
        
        TDIProfile <- cpp_TAI(as.matrix(ExpressionMatrix),as.vector(Divergencestratum))
        names(TDIProfile) <- names(ExpressionMatrix)
        
        return(TDIProfile)
        
}