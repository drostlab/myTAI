context("Test: REMatrix() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

equal_df <- function(df1,df2){
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1,df2))
}


TestREMatrix1 <- t(apply(dplyr::summarise_each(dplyr::group_by(PhyloExpressionSetExample[,c(1,3:9)],Phylostratum),dplyr::funs(mean))[ , 2:8],1,function(x) (x-min(x))/(max(x)-min(x))))
rownames(TestREMatrix1) <- 1:12

TestREMatrix2 <- t(apply(dplyr::summarise_each(dplyr::group_by(DivergenceExpressionSetExample[,c(1,3:9)],Divergence.stratum),dplyr::funs(mean))[ , 2:8],1,function(x) (x-min(x))/(max(x)-min(x))))
rownames(TestREMatrix2) <- 1:10


nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to REMatrix()",{
        expect_error(REMatrix(nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})


test_that("REMatrix() computes correct values...",{
        expect_true(equal_df(REMatrix(PhyloExpressionSetExample), TestREMatrix1))
        expect_true(equal_df(REMatrix(DivergenceExpressionSetExample), TestREMatrix2))
})

test_that("error occurs when when REMatrix() is computed on only one developmental stage ",{
        expect_error(REMatrix(PhyloExpressionSetExample[ , 1:3]))
})


