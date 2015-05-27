context("Test: pTDI() ")

data(DivergenceExpressionSetExample)


equal_df <- function(df1,df2){
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1,df2))
}


test_that("pTDI computes correct partial TDI contribution values...",{
        
        expect_true(equal_df(pTDI(DivergenceExpressionSetExample),apply(pStrata(DivergenceExpressionSetExample),2,cumsum)))
        
        expect_equal(round(as.vector(pTDI(DivergenceExpressionSetExample)[1, ]),7),c(0.2174378,0.2207644,0.2309211,0.2214881,0.2195601,0.2047938,0.1704023))
        
        expect_equal(as.vector(TDI(DivergenceExpressionSetExample)), as.vector(pTDI(DivergenceExpressionSetExample)[10, ]))
        
})




