context("Test : tf() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

equal_df <- function(df1,df2){
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1,df2))
}

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to tf()",{
        expect_error(tf(nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})

Test.tf.Func <- function(ExpressionSet, FUN){
        
        f <- match.fun(FUN)
        ncols <- ncol(ExpressionSet)
        res <- cbind(ExpressionSet[ , 1:2],apply(ExpressionSet[ , 3:ncols],2,f))
        return(res)
}

test_that("tf() computes correct transformation values ",{
        expect_true(equal_df(tf(PhyloExpressionSetExample , log2), Test.tf.Func(PhyloExpressionSetExample,log2)))
        expect_true(equal_df(tf(DivergenceExpressionSetExample , log2), Test.tf.Func(DivergenceExpressionSetExample,log2)))
        
        expect_true(equal_df(tf(PhyloExpressionSetExample , sqrt), Test.tf.Func(PhyloExpressionSetExample,sqrt)))
        expect_true(equal_df(tf(DivergenceExpressionSetExample , sqrt), Test.tf.Func(DivergenceExpressionSetExample,sqrt)))
        
        expect_true(equal_df(tf(PhyloExpressionSetExample , log2), Test.tf.Func(PhyloExpressionSetExample,log2)))
        expect_true(equal_df(tf(DivergenceExpressionSetExample , sqrt), Test.tf.Func(DivergenceExpressionSetExample,sqrt)))
        
        expect_true(equal_df(tf(PhyloExpressionSetExample , function(x) x/2), Test.tf.Func(PhyloExpressionSetExample,function(x) x/2)))
        expect_true(equal_df(tf(DivergenceExpressionSetExample , function(x) x/2), Test.tf.Func(DivergenceExpressionSetExample,function(x) x/2)))
})







