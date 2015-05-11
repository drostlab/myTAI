context("Test: FilterRNASeqCT() ")

data(PhyloExpressionSetExample)

# adapted from: https://github.com/hadley/dplyr/blob/master/tests/testthat/test-arrange.r 
equal_df <- function(df1, df2) {
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1, df2))
}

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to FilterRNASeqCT()",{
        expect_error(FilterRNASeqCT(ExpressionSet = nonStandardExpressionSet,
                                    cut.off       = 1000),"The present input object does not fulfill the ExpressionSet standard.")
})


# a test set for a complete PhyloExpressionSet in ExpressionSet notation (standard)
TestExpressionSet_completePES <- PhyloExpressionSetExample[1:10, ]

test_that("correct rows are removed (filtered) from the count table. Method: 'const'",{
        
        expect_true(equal_df(FilterRNASeqCT(TestExpressionSet_completePES,1000,"const"),TestExpressionSet_completePES[-c(1,3,4,6,8,9), ]))
        
})


test_that("correct rows are removed (filtered) from the count table. Method: 'min-set'",{
        
        expect_true(equal_df(FilterRNASeqCT(TestExpressionSet_completePES,1000,"min-set"),TestExpressionSet_completePES[-c(3,6,8,9), ]))
        
})


test_that("correct rows are removed (filtered) from the count table. Method: 'n-set'",{
        
        expect_true(equal_df(FilterRNASeqCT(TestExpressionSet_completePES,800,"n-set",5),TestExpressionSet_completePES[-c(6,8,9), ]))
        
})


test_that("error occurs when n is larger than the number of available stages when choosing method = 'n-set'",{
        
        expect_error(FilterRNASeqCT(TestExpressionSet_completePES,800,"n-set",8),"n is larger than the number of available stages in your ExpressionSet...")
})


test_that("error occurs when method = 'n-set', but n = NULL",{
        
        expect_error(FilterRNASeqCT(TestExpressionSet_completePES,800,"n-set"),"Please specify the number of stages n for which expresssion levels need to be above the cutoff to be retained in the count table.")
})



