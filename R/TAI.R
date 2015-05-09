#' @title Compute the Transcriptome Age Index (TAI)
#' @description
#' This function computes the phylogenetically based transcriptome age index (TAI) introduced by Domazet-Loso & Tautz, 2010.
#' @param PhyloExpressionSet a standard PhyloExpressionSet object.
#' @details The TAI measure represents the weighted arithmetic mean (expression levels as
#' weights for the phylostratum value) over all evolutionary age categories denoted as \emph{phylostra}.
#'
#' \deqn{TAI_s = \sum (e_is * ps_i) / \sum e_is}
#'
#' where TAI_s denotes the TAI value in developmental stage s, 
#' e_is denotes the gene expression level of gene i in stage s, 
#' and ps_i denotes the corresponding phylostratum of gene i, \eqn{i = 1,...,N} and N = total number of genes.
#'
#'
#' Internally the function calls the C++ function \code{cpp_TAI} to speed up TAI computations.
#' 
#' @return a numeric vector containing the TAI values for all given developmental stages.
#' @references 
#' Domazet-Loso T. and Tautz D. (2010). \emph{A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns}. Nature (468): 815-818.
#'
#' Quint M et al. (2012). \emph{A transcriptomic hourglass in plant embryogenesis}. Nature (490): 98-101.
#' 
#' Drost HG et al. (2015) \emph{Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis}. Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012.
#' 
#' @author Hajk-Georg Drost
#' @seealso \code{\link{TDI}}, \code{\link{PlotPattern}}, \code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}}
#' @examples
#' 
#' # reading a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # computing the TAI profile of a given PhyloExpressionSet object
#' TAIs <- TAI(PhyloExpressionSetExample)
#' 
#' 
#' 
#' @useDynLib myTAI
#' @importFrom Rcpp sourceCpp
#' @export

TAI <- function(PhyloExpressionSet)
{
        
        is.ExpressionSet(PhyloExpressionSet)
        
        nCols <- dim(PhyloExpressionSet)[2]
        ExpressionMatrix <- PhyloExpressionSet[ , 3:nCols]
        Phylostratum <- PhyloExpressionSet[ , 1]
        TAIProfile <- vector(mode = "numeric",length = nCols-2)
        
        TAIProfile <- cpp_TAI(as.matrix(ExpressionMatrix),as.vector(Phylostratum))
        names(TAIProfile) <- names(ExpressionMatrix)
        
        return(TAIProfile)
        
}