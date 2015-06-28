context("Test: Expressed() ")

data(PhyloExpressionSetExample)

# adapted from: https://github.com/hadley/dplyr/blob/master/tests/testthat/test-arrange.r 
equal_df <- function(df1, df2) {
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1, df2))
}

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to Expressed()",{
        expect_error(Expressed(ExpressionSet = nonStandardExpressionSet,
                                    cut.off       = 1000),
                     "The present input object does not fulfill the ExpressionSet standard.")
})


test_that("error occurs when wrong filter method is specified in Expressed()",{
        
        expect_error(Expressed(TestExpressionSet_completePES,1000,"m-set"),
                     "Please specify a filter method that is implemented in this function!")
        
})

# a test set for a complete PhyloExpressionSet in ExpressionSet notation (standard)
TestExpressionSet_completePES <- PhyloExpressionSetExample[1:10, ]

test_that("correct rows are removed (filtered) from the count table. Method: 'const' and comparison : 'below'",{
        
        expect_true(equal_df(Expressed(TestExpressionSet_completePES,1000,"const"),
                             TestExpressionSet_completePES[-c(1,3,4,6,8,9), ]))
        
})


test_that("correct rows are removed (filtered) from the count table. Method: 'const' and comparison : 'above'",{
        
        expect_true(equal_df(Expressed(TestExpressionSet_completePES,1000,"const",comparison = "above"),
                             TestExpressionSet_completePES[c(6,8,9), ]))
        
})


test_that("correct rows are removed (filtered) from the count table. Method: 'const' and comparison : 'both'",{
        
        expect_true(equal_df(Expressed(TestExpressionSet_completePES,c(800,2000),"const",comparison = "both"),
                             TestExpressionSet_completePES[c(2,3,4,7), ]))
        
})


test_that("error occurs when no genes fulfill the threshold criteria: comparison = 'below'",{
        
        expect_error(Expressed(TestExpressionSet_completePES,500,"const",comparison = "below"),
                     "None of the genes fulfilles the threshold criteria. Please choose a less conservative threshold or filter method.")
        
})

test_that("error occurs when no genes fulfill the threshold criteria: comparison = 'above'",{
        
        expect_error(Expressed(TestExpressionSet_completePES,80000,"const",comparison = "above"),
                     "None of the genes fulfilles the threshold criteria. Please choose a less conservative threshold or filter method.")
        
})

test_that("error occurs when no genes fulfill the threshold criteria: comparison = 'both'",{
        
        expect_error(Expressed(TestExpressionSet_completePES,c(1000,1100),"const",comparison = "both"),
                             "None of the genes fulfilles the threshold criteria. Please choose a less conservative threshold or filter method.")
        
})


test_that("correct rows are removed (filtered) from the count table. Method: 'min-set'; comparison = 'below'",{
        
        expect_true(equal_df(Expressed(TestExpressionSet_completePES,1000,"min-set",comparison = "below"),
                             TestExpressionSet_completePES[-c(3,6,8,9), ]))
        
})

test_that("correct rows are removed (filtered) from the count table. Method: 'min-set'; comparison = 'above'",{
        
        expect_true(equal_df(Expressed(TestExpressionSet_completePES,1000,"min-set",comparison = "above"),
                             TestExpressionSet_completePES[c(3,6,8,9), ]))
        
})


test_that("correct rows are removed (filtered) from the count table. Method: 'min-set'; comparison = 'both'",{
        
        expect_true(equal_df(Expressed(TestExpressionSet_completePES,c(800,2000),"min-set",comparison = "both"),
                             TestExpressionSet_completePES[c(1:4,6,7), ]))
        
})


test_that("correct rows are removed (filtered) from the count table. Method: 'n-set'; comparison = 'below'",{
        
        expect_true(equal_df(Expressed(TestExpressionSet_completePES,800,"n-set",comparison = "below",n = 5),
                             TestExpressionSet_completePES[-c(6,8,9), ]))
        
})


test_that("correct rows are removed (filtered) from the count table. Method: 'n-set'; comparison = 'above'",{
        
        expect_true(equal_df(Expressed(TestExpressionSet_completePES,2000,"n-set",comparison = "above",n = 3),
                             TestExpressionSet_completePES[c(1:4,6:9), ]))
        
})


test_that("correct rows are removed (filtered) from the count table. Method: 'n-set'; comparison = 'both'",{
        
        expect_true(equal_df(Expressed(TestExpressionSet_completePES,c(900,2000),"n-set",comparison = "both",n = 4),
                             TestExpressionSet_completePES[c(1,2,4,7), ]))
        
})

test_that("error occurs when n is larger than the number of available stages when choosing method = 'n-set'",{
        
        expect_error(Expressed(TestExpressionSet_completePES,800,"n-set",n = 8),
                     "n is larger than the number of available stages in your ExpressionSet...")
})


test_that("error occurs when method = 'n-set', but n = NULL",{
        
        expect_error(Expressed(TestExpressionSet_completePES,800,"n-set"),
                     "Please specify the number of stages n for which expresssion levels need to be above the cutoff to be retained in the count table.")
})


test_that("n-set with n = 4 computes the same values as min-set for 7 stages; comparison = 'below'",{
        
        expect_true(equal_df(Expressed(PhyloExpressionSetExample,8000,"n-set",comparison = "below",n = 4),
                             Expressed(PhyloExpressionSetExample,8000,"min-set",comparison = "below")))
        
})

test_that("n-set with n = 4 computes the same values as min-set for 7 stages; comparison = 'above'",{
        
        expect_true(equal_df(Expressed(PhyloExpressionSetExample,2000,"n-set",comparison = "above",n = 4),
                             Expressed(PhyloExpressionSetExample,2000,"min-set",comparison = "above")))
        
})

test_that("n-set with n = 4 computes the same values as min-set for 7 stages; comparison = 'both'",{
        
        expect_true(equal_df(Expressed(PhyloExpressionSetExample,c(900,2000),"n-set",comparison = "both",n = 4),
                             Expressed(PhyloExpressionSetExample,c(900,2000),"min-set",comparison = "both")))
        
})

test_that("n-set with n = 6 computes the same values as const; comparison = 'below'",{
        
        expect_true(equal_df(Expressed(PhyloExpressionSetExample,2000,"n-set",comparison = "below",n = 6),
                             Expressed(PhyloExpressionSetExample,2000,"const",comparison = "below")))
        
})

test_that("n-set with n = 6 computes the same values as const; comparison = 'above'",{
        
        expect_true(equal_df(Expressed(PhyloExpressionSetExample,8000,"n-set",comparison = "above",n = 6),
                             Expressed(PhyloExpressionSetExample,8000,"const",comparison = "above")))
        
})

test_that("n-set with n = 6 computes the same values as const; comparison = 'both'",{
        
        expect_true(equal_df(Expressed(PhyloExpressionSetExample,c(900,2000),"n-set",comparison = "both",n = 6),
                             Expressed(PhyloExpressionSetExample,c(900,2000),"const",comparison = "both")))
        
})


test_that("error occurs when only one cut.off is specified when selecting comparison = 'both'",{
        
        expect_error(Expressed(PhyloExpressionSetExample,2000,"n-set",comparison = "both",n = 6),
                     "When choosing: comparison == 'both', the cut.off argument needs to store two cut.off values: lower-cut.off and upper-cut.off")
})


