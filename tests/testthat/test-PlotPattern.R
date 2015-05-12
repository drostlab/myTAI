context("Test: PlotPattern() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotPattern()",{
        expect_error(PlotPattern(nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})
