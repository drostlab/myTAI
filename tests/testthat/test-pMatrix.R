context("Test: pMatrix() ")


data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

equal_df <- function(df1,df2){
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1,df2))
}

partialTAI <- function(ExpressionSet){
        
        ncols <- ncol(ExpressionSet)
        #s <- colSums(ExpressionSet[ , 3:ncols])
        partialTXI <- matrix(NA_real_,nrow = nrow(ExpressionSet),ncol = ncols-2)
        for(i in 1:(ncols-2)){
                partialTXI[ , i] <- ExpressionSet[ , 1] * (ExpressionSet[ , i+2]/sum(ExpressionSet[ , i+2]))
        }
        
        colnames(partialTXI) <- colnames(ExpressionSet)[3:ncols]
        rownames(partialTXI) <- ExpressionSet[ , 2]
        return(partialTXI)
}


nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to pMatrix()",{
        expect_error(pMatrix(nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})


test_that("pMatrix() computes correct partial TAI values ",{
        
        expect_true(equal_df(pMatrix(PhyloExpressionSetExample),partialTAI(PhyloExpressionSetExample)))
})

test_that("pMatrix() computes correct partial TDI values ",{
        
        expect_true(equal_df(pMatrix(DivergenceExpressionSetExample),partialTAI(DivergenceExpressionSetExample)))
})

test_that("colSums(pMatrix(ExpressionSet)) == TAI(ExpressionSet)",{
        
        expect_equal(colSums(pMatrix(PhyloExpressionSetExample)),TAI(PhyloExpressionSetExample))
        expect_equal(colSums(pMatrix(DivergenceExpressionSetExample)),TDI(DivergenceExpressionSetExample))
})


test_that("colSums(pMatrix(ExpressionSet)) == colSums(pStrata(ExpressionSet))",{
        
        expect_equal(colSums(pMatrix(PhyloExpressionSetExample)),colSums(pStrata(PhyloExpressionSetExample)))
        expect_equal(colSums(pMatrix(DivergenceExpressionSetExample)),colSums(pStrata(DivergenceExpressionSetExample)))
        
})

