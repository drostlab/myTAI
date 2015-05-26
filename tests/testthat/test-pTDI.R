context("Test: pTDI() ")

data(DivergenceExpressionSetExample)


equal_df <- function(df1,df2){
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1,df2))
}


test_that("pTDI computes correct partial TDI contribution values...",{
        
        expect_true(equal_df(pTAI(DivergenceExpressionSetExample),apply(pStrata(DivergenceExpressionSetExample),2,cumsum)))
        
        
})




