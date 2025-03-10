
# expressed genes

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