context("Test: PlotMeans() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotMeans()",{
        expect_error(PlotMeans(nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})


test_that("error occurs when Groups are not specified ",{
        expect_error(PlotMeans(PhyloExpressionSetExample),"Your Groups list does not store any items.")
})

test_that("error occurs when legendName is not specified",{
        
        expect_error(PlotMeans(PhyloExpressionSetExample,list(1:3,4:12)),"lease specify the type of ExpressionSet you are working with: legendName = 'PS' or 'DS'.")
})

test_that("error occurs when PS or DS are specified in Groups which are not present in the ExpressionSet ",{
        
        
        expect_error(PlotMeans(PhyloExpressionSetExample,list(1:3,4:13),"PS"),"There are items in your Group elements that are not available in the age column of your ExpressionSet.")
        expect_error(PlotMeans(DivergenceExpressionSetExample,list(1:12),"DS"),"There are items in your Group elements that are not available in the age column of your ExpressionSet.")
        
})

