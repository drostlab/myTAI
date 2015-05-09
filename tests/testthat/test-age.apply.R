context("Test: age.apply() ")

# adapted from: https://github.com/hadley/dplyr/blob/master/tests/testthat/test-arrange.r 
equal_df <- function(df1, df2) {
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1, df2))
}

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to age.apply()",{
        expect_error(age.apply(nonStandardExpressionSet , RE),"The present input object does not fulfill the ExpressionSet standard.")
})

test_that("age.apply() + RE equals REMatrix() for PhyloExpressionSetExample", {
        expect_true(equal_df(age.apply(PhyloExpressionSetExample, RE)
, REMatrix(PhyloExpressionSetExample)))
        
})

test_that("age.apply() + RE equals REMatrix() for DivergenceExpressionSetExample", {
        expect_true(equal_df(age.apply(DivergenceExpressionSetExample, RE)
                             , REMatrix(DivergenceExpressionSetExample)))
        
})







