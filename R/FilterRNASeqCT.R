#' @title Filter Expression Levels in RNASeq Count Tables
#' @description This function takes an ExpressionSet object that is based
#' on a RNASeq count table (CT) and removes genes from this CT that
#' have an expression level below a defined \code{cut.off} value.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param cut.off a numeric value specifying the expression cut off to define genes as \emph{not expressed} in case they express values that lie below this threshold.
#' @param method a method how to treat gene expression in multiple stages. Options are \code{"const"}, \code{"min-set"}, and \code{"n-set"}.
#' @param n a numeric value for \code{method = "n-set"}.
#' @author Hajk-Georg Drost
#' @details
#' 
#' The
#' @examples
#' 
#' # 
#' 
#' @export

FilterRNASeqCT <- function(ExpressionSet, cut.off, method = "const", n = NULL){
        
        is.ExpressionSet(ExpressionSet)
        
        colnames(ExpressionSet)[1] <- "GeneID"
        ncols <- ncol(ExpressionSet)
        
        if(method == "const"){
                # determine non expressed genes (NEGs)
                NEGs <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off))) ) )
        }
        else if (method == "min-set"){
                CandidateSet <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off))) ) )
# count for each gene how many stages are above the cutoff; aco = above cut off
MinSet <- dplyr::summarise(dplyr::group_by(ExpressionSet[CandidateSet, ], GeneID), aco = function(x) sum(which(x > cutoff)))
MinSet.Filtered <- dplyr::filter(MinSet, aco <= ceiling((ncols-2)/2))
NEGs <- as.vector(dplyr::select(MinSet.Filtered, GeneID))
        } 
else if (method == "n-set"){
        
        if(is.null(n))
                stop("Please specify the number of stages n for which expresssion levels need to be above the cutoff to be retained in the count table.")
        
        CandidateSet <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off))) ) )
# count for each gene how many stages are above the cutoff; aco = above cut off
MinSet <- dplyr::summarise(dplyr::group_by(ExpressionSet[CandidateSet, ], GeneID), aco = function(x) sum(which(x > cutoff)))
MinSet.Filtered <- dplyr::filter(MinSet, aco <= n)
NEGs <- as.vector(dplyr::select(MinSet.Filtered, GeneID))
} 

if(length(NEGs) > 1){
        return(ExpressionSet[ -NEGs , ])
        
} else {
        return(ExpressionSet)
}
}