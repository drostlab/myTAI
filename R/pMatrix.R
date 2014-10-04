# Computing the Partial TAI Values for each gene and each developmental stage
# The sum over all Partial TAI Values = the total TAI value (TAIs)
pMatrix <- function(ExpressionSet)
{
        
        is.ExpressionSet(ExpressionSet)
        
        nCols <- dim(ExpressionSet)[2]
        nRows <- dim(ExpressionSet)[1]
        pTAIMatrix <- matrix(nrow = nRows,ncol = nCols-2)
        
        pTAIMatrix <- cpp_pMatrix(as.matrix(ExpressionSet[ , 3:nCols]),as.vector(ExpressionSet[ , 1]))
        
        colnames(pTAIMatrix) <- names(ExpressionSet)[3:nCols]
        rownames(pTAIMatrix) <- ExpressionSet[ , 2]
        
        return(pTAIMatrix)
        
}

