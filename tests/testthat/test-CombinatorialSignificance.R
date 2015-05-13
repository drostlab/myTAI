context("Test: CombinatorialSignificance() ")

data(PhyloExpressionSetExample)


nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to CombinatorialSignificance()",{
        expect_error(CombinatorialSignificance(ExpressionSet = nonStandardExpressionSet,
                                               replicates    = c(2,3,2),
                                               TestStatistic = "FlatLineTest",
                                               permutations  = 10,
                                               parallel      = FALSE),"The present input object does not fulfill the ExpressionSet standard.")
})

test_that("variable and constant replicate values for corresponding stages work properly..", {
        
        expect_equal(length(CombinatorialSignificance(ExpressionSet = PhyloExpressionSetExample,
                                                      replicates    = c(2,3,2),
                                                      TestStatistic = "FlatLineTest",
                                                      permutations  = 10,
                                                      parallel      = FALSE)), prod(c(2,3,2)))
        
        expect_equal(length(CombinatorialSignificance(ExpressionSet = PhyloExpressionSetExample,
                                                      replicates    = c(3,1,3),
                                                      TestStatistic = "FlatLineTest",
                                                      permutations  = 10,
                                                      parallel      = FALSE)), prod(c(3,1,3)))
        
        expect_equal(length(CombinatorialSignificance(ExpressionSet = PhyloExpressionSetExample[ , 1:8],
                                                      replicates    = 2,
                                                      TestStatistic = "FlatLineTest",
                                                      permutations  = 10,
                                                      parallel      = FALSE)), 8)
        
#         expect_error(CombinatorialSignificance(ExpressionSet = PhyloExpressionSetExample,
#                                                replicates    = c(0,3,4),
#                                                TestStatistic = "FlatLineTest",
#                                                permutations  = 10,
#                                                parallel      = FALSE))
        
})


test_that("Error occurs in case a TestStatistic other than FlatLineTest has been chosen...", {
        
        expect_error(CombinatorialSignificance(ExpressionSet = PhyloExpressionSetExample,
                                  replicates    = c(2,3,2),
                                  TestStatistic = "ReductiveHourglassTets",
                                  permutations  = 10,
                                  parallel      = FALSE), "Please enter a correct string for the test statistic: 'FlatLineTest'.")
})


