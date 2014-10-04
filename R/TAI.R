# Transcriptome age index
TAI <- function(PhyloExpressionSet)
{
        
        is.ExpressionSet(PhyloExpressionSet)
        
        nCols <- dim(PhyloExpressionSet)[2]
        ExpressionMatrix <- PhyloExpressionSet[ , 3:nCols]
        Phylostratum <- PhyloExpressionSet[ , 1]
        TAIProfile <- vector(mode = "numeric",length = nCols-2)
        
        TAIProfile <- cpp_TAI(as.matrix(ExpressionMatrix),as.vector(Phylostratum))
        names(TAIProfile) <- names(ExpressionMatrix)
        
        return(TAIProfile)
        
}