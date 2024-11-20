#' @title Top quantile expression genes
#' @description Select genes residing in the top quantile according to the mean of their expression across the stages
#'
#' @param ExpressionSet A standard ExpressionSet
#' @param p The quantile probability. Default is \code{p = .99}
#'
#' @return a character vector containing the gene ids residing in the top expression quantile
#' @export
#'
#' @examples
#' # reading a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # select genes with highest variance (top 2%)  
#' genes.top_expression <- TopExpressionGenes(PhyloExpressionSetExample, p=.98)
#' 
#' # remove top genes from the PhyloExpressionSet
#' PhyloExpressionSet.top_removed <- subset(PhyloExpressionSetExample, 
#'                                          !(GeneID %in% genes.top_expression))
#' 
#' # plot TAI of set with removed quantile
#' PlotSignature(ExpressionSet = PhyloExpressionSet.top_removed,
#'                 measure       = "TAI", 
#'                 TestStatistic = "FlatLineTest",
#'                 xlab          = "Ontogeny", 
#'                 ylab          = "TAI" )
#'
#' @author Stefan Manolache
#' 
TopExpressionGenes <- function(ExpressionSet, p = .99){
    is.ExpressionSet(ExpressionSet)
    
    # For each gene, get the average expression level across the process
    ExprProfile <- ExpressionSet[ , -1]
    ExprProfile$AvgExpression <- rowMeans(ExprProfile[, -1])
    
    
    # Select genes which are above the specified probability
    ExprProfile.top <- ExprProfile[stats::quantile(ExprProfile$AvgExpression, 
                                        probs=c(p)) < ExprProfile$AvgExpression, ]
    
    return(ExprProfile.top$GeneID)
    
}