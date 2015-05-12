context("Test: RE() ")


data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)


testRE <- function(x){
        
        y <- colMeans(x)
        return( (y-min(y))/(max(y) - min(y)) )
}

test_that("RE() computes correct values", {
        expect_true(isTRUE(all.equal(RE(PhyloExpressionSetExample[ which(PhyloExpressionSetExample[ , 1] == 1), 3:9 ]), testRE(PhyloExpressionSetExample[ which(PhyloExpressionSetExample[ , 1] == 1), 3:9 ]))))
        
        
        expect_true(isTRUE(all.equal(RE(PhyloExpressionSetExample[ which(PhyloExpressionSetExample[ , 1] == 5), 3:9 ]), testRE(PhyloExpressionSetExample[ which(PhyloExpressionSetExample[ , 1] == 5), 3:9 ]))))
        
        
        expect_true(isTRUE(all.equal(RE(PhyloExpressionSetExample[ which(PhyloExpressionSetExample[ , 1] == 12), 3:9 ]), testRE(PhyloExpressionSetExample[ which(PhyloExpressionSetExample[ , 1] == 12), 3:9 ]))))
})



