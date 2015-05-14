#' @title Transform to Relative Expression Levels
#' @description This function computes the relative expression profiles of any given gene expression set. 
#' The relative expression profile is being computed as follows:
#' 
#' \deqn{f_s = ( e_s - e_min ) / ( e_max - e_min )}
#'
#' where \eqn{e_min} and \eqn{e_max} denote the minimum/maximum mean expression level
#' over the  developmental stages s. This linear transformation corresponds to a
#' shift by \eqn{e_min} and a subsequent shrinkage by \eqn{e_max - e_min}. 
#' As a result, the relative expression level \eqn{f_s} of developmental stage s
#' with minimum \eqn{e_s} is 0, the relative expression level \eqn{f_s} of the
#' developmental stage s with maximum \eqn{e_s} is 1, and the relative expression
#' levels \eqn{f_s} of all other stages s range between 0 and 1, accordingly.
#' @param ExpressionMatrix a numeric matrix representing a gene expression matrix for which the relative expression profile shall be computed.
#' @return a vector containing the relative expression profile of the correspnding data matrix.
#' @references 
#' Domazet-Loso T and Tautz D. (2010). \emph{A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns}. Nature (468): 815-818.
#'
#' Quint M et al. (2012). \emph{A transcriptomic hourglass in plant embryogenesis}. Nature (490): 98-101.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{REMatrix}}, \code{\link{PlotRE}}
#' @examples
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#'
#' # relative expression profile of PS1 genes
#' RE(PhyloExpressionSetExample[ which(PhyloExpressionSetExample[ , 1] == 1), 3:9 ])
#'
#' 
#' @export

RE <- function(ExpressionMatrix)
{
        mDimensions <- dim(ExpressionMatrix)
        cMeans <- vector(mode = "numeric",length=mDimensions[2])
        cMeans <- colMeans(ExpressionMatrix)
        
        f_max <- max(cMeans)
        f_min <- min(cMeans)
        RE <- (cMeans - f_min) / (f_max - f_min)
        return(RE)
}

#' @title Compute a Relative Expression Matrix
#' @description This function computes the relative expression profiles of 
#' all given phylostrata or divergence-strata within a given PhyloExpressionSet or DivergenceExpressionSet.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @details For each phylostratum or divergence-stratum the corresponding relative expression profile is being computed as follows:
#'
#' \deqn{f_js = ( e_js - e_j min ) / ( e_j max - e_j min )}
#'
#' where \eqn{e_j min} and \eqn{e_j max} denote the minimum/maximum mean expression level of 
#' phylostratum j over the  developmental stages s. This linear transformation corresponds
#' to a shift by \eqn{e_j min} and a subsequent shrinkage by \eqn{e_j max - e_j min}. 
#' As a result, the relative expression level \eqn{f_js} of developmental stage s with minimum \eqn{e_js} is 0,
#' the relative expression level \eqn{f_js} of the developmental stage s with maximum \eqn{e_js} is 1, 
#' and the relative expression levels \eqn{f_js} of all other stages s range between 0 and 1, accordingly.
#' @author Hajk-Georg Drost
#' @references
#' Domazet-Loso T and Tautz D. (2010). \emph{A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns}. Nature (468): 815-818.
#'
#' Quint M et al. (2012). \emph{A transcriptomic hourglass in plant embryogenesis}. Nature (490): 98-101.
#' 
#' @seealso \code{\link{RE}}, \code{\link{PlotRE}}, \code{\link{PlotBarRE}}
#' @examples
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#'
#' # example PhyloExpressionSet
#' REMatrix(PhyloExpressionSetExample)
#'
#' # example DivergenceExpressionSet
#' REMatrix(DivergenceExpressionSetExample)
#' 
#' 
#' @export

REMatrix <- function(ExpressionSet)
{
        if(ncol(ExpressionSet) < 4)
                stop("You need at least 2 stages (= expression columns in your data.frame) to comptute relative expression levels...")
        is.ExpressionSet(ExpressionSet)
        return(age.apply(ExpressionSet = ExpressionSet, RE))
        
}

