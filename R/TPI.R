#' @title Compute the Transcriptome Polymorphism Index (TPI)
#' @description This function computes the Transcriptome Polymorphism Index (TPI) introduced by
#' Gossmann et al., 2015.
#' @param PolymorphismExpressionSet a standard PolymorphismExpressionSet object.
#' @details
#' 
#' The TPI measure represents the weighted arithmetic mean (expression levels as
#' weights) for the synonymous vs non-synonymous polymorphism ratios.
#'
#'
#' \deqn{TPI_s = \sum (e_is * P_N/N / ((P_S + 1) / S)) / \sum e_is}
#' 
#' where TPI_s denotes the TPI value in developmental stage s, e_is denotes the gene expression level of gene i in stage s, n denotes the number of genes, PN and PS denote the numbers of non-synonymous and synonymous polymorphisms, and N and S are the numbers of nonsynonymous and synonymous sites, respectively. 
#' 
#' Internally the function is written in C++ to speed up TPI computations.
#' @return a numeric vector containing the TPI values for all given developmental stages.
#' @references 
#' Gossmann et al. (2015). \emph{Transcriptomes of Plant Gametophytes Have a Higher Proportion of Rapidly Evolving and Young Genes than Sporophytes}. Mol Biol Evol. 33 (7): 1669-1678.
#' 
#' @author Hajk-Georg Drost
#' @seealso \code{\link{TAI}}, \code{\link{TDI}},  \code{\link{PlotSignature}}, \code{\link{PlotPattern}}, \code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}}
#' @examples
#' \dontrun{
#' # reading a standard PolymorphismExpressionSet
#' data(PolymorphismExpressionSetExample)
#'
#' # computing the TPI profile of a given PolymorphismExpressionSet object
#' TPIs <- TPI(PolymorphismExpressionSet)
#' }
#' 
#'
#' @export

TPI <- function(PolymorphismExpressionSet)
{
        
        is.ExpressionSet(PolymorphismExpressionSet)
        
        nCols <- dim(PolymorphismExpressionSet)[2]
        ExpressionMatrix <- dplyr::select(PolymorphismExpressionSet, 3:ncol(PolymorphismExpressionSet))
        Polymorphisms <- unlist(dplyr::select(PolymorphismExpressionSet, 1))
        TPIProfile <- vector(mode = "numeric", length = nCols - 2)
        
        
        TPIProfile <- cpp_TAI(as.matrix(ExpressionMatrix), as.vector(Polymorphisms))
        names(TPIProfile) <- names(ExpressionMatrix)
        
        return(TPIProfile)
        
}
