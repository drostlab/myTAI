RE <- function(ExpressionMatrix)
{
        mDimensions <- dim(ExpressionMatrix)
        cMeans <- vector(mode = "numeric",length=mDimensions[2])
        cMeans <- colMeans(ExpressionMatrix)
        
        f_max <- max(cMeans)
        f_min <- min(cMeans)
        RE <- (cMeans - f_min) / (f_max - f_min)
        return(RE)
}


REMatrix <- function(ExpressionSet)
{
        
        is.ExpressionSet(ExpressionSet)
        return(age.apply(ExpressionSet = ExpressionSet, RE))
        
}

