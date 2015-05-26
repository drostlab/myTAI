context("Test: pTAI() ")

data(PhyloExpressionSetExample)


equal_df <- function(df1,df2){
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1,df2))
}


test_that("pTAI computes correct partial TAI contribution values...",{
        
        expect_true(equal_df(pTAI(PhyloExpressionSetExample),apply(pStrata(PhyloExpressionSetExample),2,cumsum)))
        
        
})






