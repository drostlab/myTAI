context("Test: TDI() ")

data(DivergenceExpressionSetExample)


nonStandardExpressionSet <- DivergenceExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to TDI()",{
        expect_error(TDI(nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})


testTAIFunc <- function(ExpressionSet){
        ncols <- ncol(ExpressionSet)
        TAIVals <- vector("numeric",ncols-2) 
        
        for(i in 1:(ncols-2)){
                TAIVals[i] <-  sum((ExpressionSet[ , 1] * ExpressionSet[ , i+2])/sum(ExpressionSet[ , i+2]))
        }
        
        names(TAIVals) <- names(ExpressionSet)[3:ncols]
        
        return(TAIVals)
}


test_that("TDI() computes correct values ",{
        expect_equal(TDI(DivergenceExpressionSetExample),testTAIFunc(DivergenceExpressionSetExample))
})





