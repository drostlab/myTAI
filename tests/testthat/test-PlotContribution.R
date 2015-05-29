context("Test: PlotContribution() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotContribution()",{
        expect_error(PlotContribution(nonStandardExpressionSet, legendName = "PS"),"The present input object does not fulfill the ExpressionSet standard.")
})


test_that("Error occurs when legendName is not specified...",{
        
        expect_error(PlotContribution(PhyloExpressionSetExample),"Please specify whether your input ExpressionSet stores 'PS' or 'DS'.")
})









