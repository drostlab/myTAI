context("Test: pTAI() ")

data(PhyloExpressionSetExample)


equal_df <- function(df1,df2){
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1,df2))
}


test_that("pTAI computes correct partial TAI contribution values...",{
        
        expect_true(equal_df(pTAI(PhyloExpressionSetExample),apply(pStrata(PhyloExpressionSetExample),2,cumsum)))
        
        expect_equal(round(as.vector(pTAI(PhyloExpressionSetExample)[1, ]),7),c(0.3929533,0.3935308, 0.4142106, 0.4115399, 0.4216806, 0.4178302, 0.3883815))
        
        expect_equal(as.vector(TAI(PhyloExpressionSetExample)), as.vector(pTAI(PhyloExpressionSetExample)[12, ]))
        
})






