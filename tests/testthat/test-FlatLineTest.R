context("Test: FlatLineTest() ")

data(PhyloExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to FlatLineTest()",{
        expect_error(FlatLineTest(nonStandardExpressionSet,
                                           permutations = 1000),"The present input object does not fulfill the ExpressionSet standard.")
})


TestBootMatrix <- bootMatrix(PhyloExpressionSetExample, 1000)

res <- FlatLineTest(PhyloExpressionSetExample,
                    custom.perm.matrix = TestBootMatrix)

estimates <- fitdistrplus::fitdist(apply(TestBootMatrix,1,var),distr = "gamma",method = "mme") 

real_score <- var(TAI(PhyloExpressionSetExample)) 


test_that("FlatLineTest() computes correct p.values...", {
                
        expect_equal(FlatLineTest(PhyloExpressionSetExample,permutations = 1000,custom.perm.matrix = TestBootMatrix)$p.value, pgamma(real_score,shape = estimates$estimate[1],rate = estimates$estimate[2],lower.tail = FALSE))
})


test_that("FlatLineTest() computes correct std.dev...", {
        
        expect_equal(sum(FlatLineTest(PhyloExpressionSetExample,permutations = 1000,custom.perm.matrix = TestBootMatrix)$std.dev), sum(apply(TestBootMatrix,2,sd)))
})











