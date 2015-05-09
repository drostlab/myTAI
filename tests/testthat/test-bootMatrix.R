context("Test: bootMatrix() ")


data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)


nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to bootMatrix()",{
        expect_error(bootMatrix(nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})


test_that("# permutations matches nrow() of bootMatrix output...", {
        expect_equal(nrow(bootMatrix(PhyloExpressionSetExample, permutations = 1)), 1)
        expect_equal(nrow(bootMatrix(PhyloExpressionSetExample, permutations = 10)), 10)
        expect_equal(nrow(bootMatrix(PhyloExpressionSetExample, permutations = 100)), 100)
        expect_equal(nrow(bootMatrix(PhyloExpressionSetExample, permutations = 1000)), 1000)
        
})


test_that("bootMatrix() returns numeric values", {
        expect_true(is.numeric(bootMatrix(PhyloExpressionSetExample, permutations = 10)))
})



