context("Test: PlotBarRE() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotBarRE()",{
        expect_error(PlotBarRE(nonStandardExpressionSet,list(1:3,4:12)),"The present input object does not fulfill the ExpressionSet standard.")
})


test_that("error occurs when Groups are not specified ",{
        expect_error(PlotBarRE(PhyloExpressionSetExample),"Your Groups list does not store any items.")
})


test_that("error occurs when PS or DS are specified in Groups which are not present in the ExpressionSet ",{
        
        
        expect_error(PlotBarRE(PhyloExpressionSetExample,list(1:3,4:13)),"There are items in your Group elements that are not available in the age column of your ExpressionSet.")
        expect_error(PlotBarRE(DivergenceExpressionSetExample,list(1:12)),"There are items in your Group elements that are not available in the age column of your ExpressionSet.")
        
})

test_that("error occurs when Group elements store less than 2 items ",{
        
        expect_error(PlotBarRE(PhyloExpressionSetExample,list(1,2)),"Each Group class needs to store at least two items.")
})



