context("Test: omitMatrix() ")

data(PhyloExpressionSetExample)

# adapted from: https://github.com/hadley/dplyr/blob/master/tests/testthat/test-arrange.r 
equal_df <- function(df1, df2) {
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        colnames(df1) <- NULL
        colnames(df2) <- NULL
        isTRUE(all.equal(df1, df2))
}

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to omitMatrix()",{
        expect_error(omitMatrix(nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})


TestOmitMatrix <- omitMatrix(PhyloExpressionSetExample)[1:50 , ]
TestOmitMatrix2 <- matrix(NA_real_,50,7)
colnames(TestOmitMatrix2) <- names(PhyloExpressionSetExample)[3:9]
rownames(TestOmitMatrix2) <- paste0("(-) ",PhyloExpressionSetExample[1:50,2])

for(i in 1:50) 
        TestOmitMatrix2[i, ] <- TAI(PhyloExpressionSetExample[-c(i), ])

test_that("omitMatrix() computes correct TAI/TDI values",{
        
        expect_true(equal_df(TestOmitMatrix,TestOmitMatrix2))
        
})


