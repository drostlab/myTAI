# Transcriptome divergence index
TDI <- function(DivergenceExpressionSet)
{
        
        is.ExpressionSet(DivergenceExpressionSet)
        
        nCols <- dim(DivergenceExpressionSet)[2]
        ExpressionMatrix <- DivergenceExpressionSet[ , 3:nCols]
        Divergencestratum <- DivergenceExpressionSet[ , 1]
        TDIProfile <- vector(mode = "numeric",length = nCols-2)
        
        
        TDIProfile <- cpp_TAI(as.matrix(ExpressionMatrix),as.vector(Divergencestratum))
        names(TDIProfile) <- names(ExpressionMatrix)
        
        return(TDIProfile)
        
}