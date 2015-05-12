context("Test: pStrata() ")


data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

equal_df <- function(df1,df2){
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1,df2))
}

partialStrata <- function(ExpressionSet){
        
        ncols <- ncol(ExpressionSet)
        #s <- colSums(ExpressionSet[ , 3:ncols])
        partialTXI <- matrix(NA_real_,nrow = nrow(ExpressionSet),ncol = ncols-2)
        for(i in 1:(ncols-2)){
                partialTXI[ , i] <- ExpressionSet[ , 1] * (ExpressionSet[ , i+2]/sum(ExpressionSet[ , i+2]))
        }
        
        partialStrata <- matrix(NA_real_,ncol=ncols-2,nrow=length(names(table(ExpressionSet[ , 1]))))
        for(j in 1:length(names(table(ExpressionSet[ , 1])))){
               partialStrata[j, ] <- colSums(partialTXI[which(ExpressionSet[ , 1] == j) , ])
        }
        
        colnames(partialStrata) <- colnames(ExpressionSet)[3:ncols]
        rownames(partialStrata) <- names(table(ExpressionSet[ , 1]))
        return(partialStrata)
}

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to pStrata()",{
        expect_error(pStrata(nonStandardExpressionSet),"The present input object does not fulfill the ExpressionSet standard.")
})


test_that("pStrata() computes correct partial TAI values ",{
        
        expect_true(equal_df(pStrata(PhyloExpressionSetExample),partialStrata(PhyloExpressionSetExample)))
})

test_that("pStrata() computes correct partial TDI values ",{
        
        expect_true(equal_df(pStrata(DivergenceExpressionSetExample),partialStrata(DivergenceExpressionSetExample)))
})

test_that("colSums(pStrata(ExpressionSet)) == TAI(ExpressionSet)",{
        
        expect_equal(colSums(pStrata(PhyloExpressionSetExample)),TAI(PhyloExpressionSetExample))
        expect_equal(colSums(pStrata(DivergenceExpressionSetExample)),TDI(DivergenceExpressionSetExample))
})


test_that("colSums(pMatrix(ExpressionSet)) == colSums(pStrata(ExpressionSet))",{
        
        expect_equal(colSums(pMatrix(PhyloExpressionSetExample)),colSums(pStrata(PhyloExpressionSetExample)))
        expect_equal(colSums(pMatrix(DivergenceExpressionSetExample)),colSums(pStrata(DivergenceExpressionSetExample)))
        
})

