tf <- function(ExpressionSet, FUN){
        
        is.ExpressionSet(ExpressionSet)
        
        ncols <- dim(ExpressionSet)[2]
        f <- match.fun(FUN)
        
        return(data.frame(ExpressionSet[ , 1:2] , apply(ExpressionSet[ , 3:ncols] , 2 , f)))
        
        
}