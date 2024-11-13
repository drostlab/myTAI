
.RowVars <- function(x, ...) {
    return(rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2]-1))
}

#' @title Top quantile expression-variance genes
#' @description Select genes residing in the top quantile according to the variance of their expression across the stages
#'
#' @param ExpressionSet A standard ExpressionSet
#' @param p The quantile probability. Default is \code{p = .99}
#'
#' @return a character vector containing the gene ids residing in the top variance expression quantile
#' @export
#'
#' @examples
#' # reading a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # select genes with highest variance (top 5%)  
#' genes.top_variance <- TopVarianceGenes(PhyloExpressionSetExample, p=.95)
#' 
TopVarianceGenes <- function(ExpressionSet, p = .99){
    is.ExpressionSet(ExpressionSet)
    
    # For each gene, get the expression variance across the process
    ExprProfile <- ExpressionSet[ , -1]
    ExprProfile$AvgExpression <- .RowVars(ExprProfile[, -1])
    
    
    # Select genes which are above the specified probability
    ExprProfile.top <- ExprProfile[stats::quantile(ExprProfile$AvgExpression, 
                                        probs=c(p)) < ExprProfile$AvgExpression, ]
    
    return(ExprProfile.top$GeneID)
    
}