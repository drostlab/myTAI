context("Test: PlotCorrelation() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotCorrelation()",{
        expect_error(PlotCorrelation(nonStandardExpressionSet,nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})
