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
#' This filter function allows users to remove genes from the \code{ExpressionSet} object that undercut a certain expression level \code{cut.off}.
#' 
#' Following extraction criteria are implemented in this function: 
#' \itemize{
#' \item \code{const}: all genes that have at least one stage that undercuts the expression \code{cut.off} will be excluded from the \code{ExpressionSet}.
#' \item \code{min-set}: genes that undercut the expression level \code{cut.off} in \code{ceiling(n/2)} stages will be excluded from the \code{ExpressionSet}.
#' \item \code{n-set}: genes that undercut the expression level \code{cut.off} in \code{n} stages will be excluded from the \code{ExpressionSet}. 
#' }
#' 
#' @examples
#' data(PhyloExpressionSetExample)
#' 
#' # remove genes that have an expression level below 8000 
#' # in at least one developmental stage
#' FilterConst <- FilterRNASeqCT(ExpressionSet = PhyloExpressionSetExample, 
#'                               cut.off       = 8000, 
#'                               method        = "const")
#'                               
#' dim(FilterConst) # check number of retained genes
#' 
#' # remove genes that have an expression level below 8000 
#' # in at least ceiling(n/2) developmental stages (in this case: ceiling(7/2) = 4 stages)
#' FilterMinSet <- FilterRNASeqCT(ExpressionSet = PhyloExpressionSetExample, 
#'                                cut.off       = 8000, 
#'                                method        = "min-set")
#'                                
#' dim(FilterMinSet) # check number of retained genes
#' 
#' # remove genes that have an expression level below 8000 
#' # in at least n developmental stages (in this case: 5 stages)
#' FilterNSet <- FilterRNASeqCT(ExpressionSet = PhyloExpressionSetExample, 
#'                              cut.off       = 8000, 
#'                              method        = "n-set",
#'                              n             = 5)
#'                                
#' dim(FilterMinSet) # check number of retained genes
#' 
#' @export

FilterRNASeqCT <- function(ExpressionSet, cut.off, method = "const", n = NULL){
        
        is.ExpressionSet(ExpressionSet)
        
        ncols <- ncol(ExpressionSet)
        
        if(method == "const"){
                # determine non expressed genes (NEGs)
                NEGs <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off))) ) )
        }
        else if (method == "min-set"){
                CandidateSet <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off))) ) )
                # count for each gene how many stages are above the cutoff; aco = above cut off
                CandidateExpressionSet <- ExpressionSet[CandidateSet, ]
                MinSet <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x > cut.off))
                NEGs <- match(CandidateExpressionSet[which(MinSet <= ceiling((ncols-2)/2)),2],ExpressionSet[ , 2])
                
        } 
        else if (method == "n-set"){
                if(is.null(n))
                        stop("Please specify the number of stages n for which expresssion levels need to be above the cutoff to be retained in the count table.")
                
                if(n > (ncols-2))
                        stop("n is larger than the number of available stages in your ExpressionSet...")
                        
                CandidateSet <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off))) ) )
                
                CandidateExpressionSet <- ExpressionSet[CandidateSet, ]
                MinSet <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x > cut.off))
                NEGs <- match(CandidateExpressionSet[which(MinSet <= n),2],ExpressionSet[ , 2])
        } 
        
        if(length(NEGs) > 1){
                return(ExpressionSet[ -NEGs , ])
        } else {
                return(ExpressionSet)
        }
}




