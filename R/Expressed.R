#' @title Filter for Expressed Genes
#' @description This function takes an ExpressionSet object and removes genes from the gene expression matrix that
#' have an expression level below, above, or below AND above a defined \code{cut.off} value. Hence, this function allows to remove
#' genes that have been defined as \emph{not expressed} or \emph{outliers} and returns an \code{ExpressionSet} retaining only expressed genes.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param cut.off a numeric value specifying the expression cut-off to define genes as \emph{not expressed} (\code{comparison = "below"}) , \emph{outliers} (\code{comparison = "above"}), or both (\code{comparison = "both"}). See \code{comparison} for details. In case \code{comparison = "both"}, the \code{cut.off} argument must be a two dimensional vector defining the lower \code{cut.off} value at the first position and the upper \code{cut.off} value at the second position.
#' @param method a method defining how to treat gene expression values in multiple stages. The corresponding method that is chosen allows to control the stage-wise fulfillment of the threshold criteria. Options are \code{"const"}, \code{"min-set"}, and \code{"n-set"}.
#' @param comparison a character string specifying whether genes having expression levels
#'  below, above, or below AND above (both) the \code{cut.off} value should be excluded from the dataset.
#'  In case \code{comparison = "both"} is chosen, the \code{cut.off} argument must be a two dimensional vector defining the lower \code{cut.off} value at the first position and the upper \code{cut.off} value
#' at the second position. 
#' @param n a numeric value for \code{method = "n-set"}.
#' @author Hajk-Georg Drost
#' @details
#' This filter function allows users to remove genes from the \code{ExpressionSet} object that undercut or exceed a certain expression level \code{cut.off}.
#' 
#' Following extraction criteria are implemented in this function: 
#' \itemize{
#' \item \code{const}: all genes that have at least one stage that undercuts or exceeds the expression \code{cut.off} will be excluded from the \code{ExpressionSet}. Hence, for a 7 stage \code{ExpressionSet} genes passing the expression level \code{cut.off} in 6 stages will be retained in the \code{ExpressionSet}.
#' \item \code{min-set}: genes passing the expression level \code{cut.off} in \code{ceiling(n/2)} stages will be retained in the \code{ExpressionSet}, where \emph{n} is the number of stages in the \code{ExpressionSet}.
#' \item \code{n-set}: genes passing the expression level \code{cut.off} in \code{n} stages will be retained in the \code{ExpressionSet}. Here, the argument \code{n} needs to be specified.
#' }
#' 
#' @examples
#' data(PhyloExpressionSetExample)
#' 
#' # remove genes that have an expression level below 8000 
#' # in at least one developmental stage
#' FilterConst <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
#'                          cut.off       = 8000, 
#'                          method        = "const",
#'                          comparison    = "below")
#'                               
#' dim(FilterConst) # check number of retained genes
#' 
#' # remove genes that have an expression level below 8000 
#' # in at least 3 developmental stages 
#' # (in this case: ceiling(7/2) = 4 stages fulfilling the cut-off criteria)
#' FilterMinSet <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
#'                           cut.off       = 8000, 
#'                           method        = "min-set",
#'                           comparison    = "below")
#'                                
#' dim(FilterMinSet) # check number of retained genes
#' 
#' # remove genes that have an expression level below 8000 
#' # in at least 5 developmental stages (in this case: n = 2 stages fulfilling the criteria)
#' FilterNSet <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
#'                         cut.off       = 8000, 
#'                         method        = "n-set",
#'                         comparison    = "below",
#'                         n             = 2)
#'                                
#' dim(FilterMinSet) # check number of retained genes
#' 
#' 
#' 
#' # remove expression levels that exceed the cut.off criteria
#' FilterMinSet <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
#'                           cut.off       = 12000, 
#'                           method        = "min-set",
#'                           comparison    = "above")
#'                                
#' dim(FilterMinSet) # check number of retained genes
#' 
#' 
#' # remove expression levels that undercut AND exceed the cut.off criteria
#' FilterMinSet <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
#'                           cut.off       = c(8000,12000), 
#'                           method        = "min-set",
#'                           comparison    = "both")
#'                                
#' dim(FilterMinSet) # check number of retained genes
#' 
#' 
#' @export

Expressed <- function(ExpressionSet,
                      cut.off,
                      method     = "const",
                      comparison = "below",
                      n          = NULL){
        
        if(!is.element(method,c("const","min-set","n-set")))
                stop("Please specify a filter method that is implemented in this function!", call. = FALSE)
        
        if (!is.element(comparison,c("below","above","both")))
                stop("Please select an appropriate comparison method implemented in this function.", call. = FALSE)
        
        if ((length(cut.off) != 2) & (comparison == "both"))
                stop("When choosing: comparison == 'both', the cut.off argument needs to store two cut.off values: lower-cut.off and upper-cut.off", call. = FALSE)
        
        if ((length(cut.off) > 1) & (comparison != "both"))
                stop("When choosing: comparison == 'below' or 'above', the cut.off argument needs to store only one cut.off value.", call. = FALSE)
        
        is.ExpressionSet(ExpressionSet)
        ncols <- ncol(ExpressionSet)
        
        if (method == "const"){
                # determine non expressed genes (NEGs)
                
                if (comparison == "below")
                        NEGs <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off))) ) )
                
                if (comparison == "above")
                        NEGs <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x > cut.off))) ) )
                
                if (comparison == "both"){
                        
                        NEG.below <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off[1]))) ) )
                        
                        NEG.above <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x > cut.off[2]))) ) )
                        
                        NEGs <- c(NEG.below,NEG.above)
                }
                
        }
        else if (method == "min-set"){
                
                if (comparison == "below")
                        CandidateSet <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off))) ) )
                
                if (comparison == "above")
                        CandidateSet <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x > cut.off))) ) )
                
                if (comparison == "both"){
                        
                        CandidateSet.below <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off[1]))) ) )
                        
                        CandidateSet.above <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x > cut.off[2]))) ) )
                        
                        CandidateSet <- c(CandidateSet.below,CandidateSet.above)
                }
                
                if(length(CandidateSet) == 0)
                        stop("None of the genes fulfilles the threshold criteria. Please choose a less conservative threshold or filter method.", call. = FALSE)
                
                # count for each gene how many stages are above the cutoff; aco = above cut off
                CandidateExpressionSet <- ExpressionSet[CandidateSet, ]
                
                if (comparison == "below"){
                        MinSet <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x > cut.off))
                        MinSetGenes <- which(MinSet <= ceiling((ncols-2)/2))
                }
                        
                if (comparison == "above"){
                        
                        MinSet <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x < cut.off))
                        MinSetGenes <- which(MinSet <= ceiling((ncols-2)/2))
                }
                        
                if (comparison == "both"){
                        
                        MinSet.below <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x > cut.off[1]))
                        MinSetGenes.below <- which(MinSet.below <= ceiling((ncols-2)/2))
                        
                        MinSet.above <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x < cut.off[2]))
                        MinSetGenes.above <- which(MinSet.above <= ceiling((ncols-2)/2))
                        
                        MinSetGenes <- c(MinSetGenes.below,MinSetGenes.above)
                }
                        
                if(length(MinSetGenes) == 0)
                        stop("None of the genes fulfilles the threshold criteria. Please choose a less conservative threshold or filter method.", call. = FALSE)
                
                NEGs <- match(CandidateExpressionSet[MinSetGenes, 2],ExpressionSet[ , 2])
                
        } 
        else if (method == "n-set"){
                if(is.null(n))
                        stop("Please specify the number of stages n for which expresssion levels need to be above the cutoff to be retained in the count table.", call. = FALSE)
                
                if(n > (ncols-2))
                        stop("n is larger than the number of available stages in your ExpressionSet...", call. = FALSE)
                        
                if (comparison == "below")
                        CandidateSet <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off))) ) )
                
                if (comparison == "above")
                        CandidateSet <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x > cut.off))) ) )
                
                if (comparison == "both"){
                        
                        CandidateSet.below <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x < cut.off[1]))) ) )
                        
                        CandidateSet.above <- unique( unlist( apply( ExpressionSet[ , 3:ncols], 2 ,function(x) list(which(x > cut.off[2]))) ) )
                        
                        CandidateSet <- c(CandidateSet.below,CandidateSet.above)
                }
                
                if(length(CandidateSet) == 0)
                        stop("None of the genes fulfilles the threshold criteria. Please choose a less conservative threshold or filter method.", call. = FALSE)
                
                CandidateExpressionSet <- ExpressionSet[CandidateSet, ]
                
                if (comparison == "below"){
                        MinSet <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x > cut.off))
                        MinSetGenes <- which(MinSet <= n)
                }
                
                if (comparison == "above"){
                        
                        MinSet <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x < cut.off))
                        MinSetGenes <- which(MinSet <= n)
                }
                
                if (comparison == "both"){
                        
                        MinSet.below <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x > cut.off[1]))
                        MinSetGenes.below <- which(MinSet.below <= n)
                        
                        MinSet.above <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x < cut.off[2]))
                        MinSetGenes.above <- which(MinSet.above <= n)
                        
                        MinSetGenes <- c(MinSetGenes.below,MinSetGenes.above)
                        
                }
#                 MinSet <- apply(CandidateExpressionSet[ , 3:ncols],1,function(x) sum(x > cut.off))
#                 MinSetGenes <- which(MinSet <= n)
                
                if(length(MinSetGenes) == 0)
                        stop("None of the genes fulfilles the threshold criteria. Please choose a less conservative threshold or filter method.", call. = FALSE)
                
                NEGs <- match(CandidateExpressionSet[MinSetGenes , 2],ExpressionSet[ , 2])
        } 
        
        if (nrow(ExpressionSet[ -NEGs , ]) > 0){
                return(ExpressionSet[ -NEGs , ])
        } else {
                stop("None of the genes fulfilles the threshold criteria. Please choose a less conservative threshold or filter method.", call. = FALSE)        
                }
}




