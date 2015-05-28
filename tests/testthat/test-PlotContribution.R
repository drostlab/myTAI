context("Test: PlotContribution() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotContribution()",{
        expect_error(PlotContribution(nonStandardExpressionSet, legendName = "PS"),"The present input object does not fulfill the ExpressionSet standard.")
})