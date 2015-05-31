context("Test: PlotEnrichment() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

set.seed(123)
test_set <- sample(PhyloExpressionSetExample[ , 2],22000)


nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotEnrichment()",{
        expect_error(PlotEnrichment(nonStandardExpressionSet,test_set,legendName = "PS"),"The present input object does not fulfill the ExpressionSet standard.")
})

test_that("error pccurs when test.set inlcudes more genes than are available in the
          input ExpressionSet",{
                  
                  expect_error(PlotEnrichment(PhyloExpressionSetExample[sample(1:22000,19000), ],test_set,legendName = "PS") , "Your input GeneID vector stores more elements than are available in your ExpressionSet object...")
                  
          })



