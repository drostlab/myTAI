#' @title Compute the Phylostratum Contribution to the Global Transcriptome Age Index
#' @description This function takes a standard \emph{ExpressionSet} object and computes the partial contribution of the different phylostrata (ps) to the global Transcriptome Age Index profile.
#' @param ExpressionSet a standard PhyloExpressionSet object.
#' @author Hajk-Georg Drost
#' @details This way of computing the partial contribution of the different phylostrata (ps) to the global Transcriptome Age Index profile was introduced by Domazet-Loso and Tautz, 2010. This function (\code{pTAI}) computes the partial TAI contribution for each phylostratum and each developmental stage and returns a data matrix storing the partial TAI contribution value for each phylostratum and each developmental stage. 
#' 
#' @examples
#' 
#' data(PhyloExpressionSetExample)
#' 
#' # get the partial contribution of phylostrata to the global
#' # TAI pattern
#' pTAI(PhyloExpressionSetExample)
#' 
#' @references 
#' Domazet-Loso T. and Tautz D. (2010). \emph{A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns}. Nature (468): 815-818.
#' @seealso \code{\link{pTDI}}, \code{\link{TAI}}, \code{\link{TDI}}, \code{\link{PlotPattern}}
#' @export
pTAI <- function(ExpressionSet){
        
        apply(pStrata(ExpressionSet),2,cumsum)
        
}




