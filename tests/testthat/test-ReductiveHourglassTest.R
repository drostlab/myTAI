context("Test: ReductiveHourglassTest() ")

data(PhyloExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to ReductiveHourglassTest()",{
        expect_error(ReductiveHourglassTest(nonStandardExpressionSet,
                                           modules = list(early = 1:2, mid = 3:5, late = 6:7),
                                           permutations = 1000),"The present input object does not fulfill the ExpressionSet standard.")
})


test_that("lillie.test is NA..", {
        
        expect_true(is.na(ReductiveHourglassTest(PhyloExpressionSetExample,
                                                modules = list(early = 1:2, mid = 3:5, late = 6:7),
                                                permutations = 1000)$lillie.test))
})

test_that("lillie.test is computed...",{
        
        expect_true(ReductiveHourglassTest(PhyloExpressionSetExample,
                                          modules = list(early = 1:2, mid = 3:5, late = 6:7),
                                          permutations = 1000, lillie.test = TRUE)$lillie.test %in% c(TRUE,FALSE))
        
})

test_that("error occurs when module selection does not match number of developmental stages..",{
        
        expect_error(ReductiveHourglassTest(PhyloExpressionSetExample,
                                           modules = list(early = 1:2, mid = 3:5, late = 6:8),
                                           permutations = 1000),"The number of stages classified into the three modules does not match the total number of stages stored in the given ExpressionSet.")
})


test_that("error occurs when modules aren't specified...",{
        
        expect_error(ReductiveHourglassTest(PhyloExpressionSetExample,
                                           permutations = 1000))
})


TestBootMatrix <- bootMatrix(PhyloExpressionSetExample, 1000)

res <- ReductiveHourglassTest(PhyloExpressionSetExample,
                             modules = list(early = 1:2, mid = 3:5, late = 6:7),
                             custom.perm.matrix = TestBootMatrix)

estimates <- fitdistrplus::fitdist(apply(TestBootMatrix,1,rhScore,early = 1:2,mid=3:5,late=6:7,method = "min",scoringMethod = "mean-mean"),distr = "norm", method = "mme") 

real_score <- rhScore(TAI(PhyloExpressionSetExample),1:2,3:5,6:7,
                      method = "min",scoringMethod = "mean-mean") 

test_that("ReductiveHourglassTest() computes correct std.dev and p.values values...", {
        
        expect_equal(res$p.value,pnorm(real_score,mean = estimates$estimate[1],sd = estimates$estimate[2],lower.tail = FALSE))
        expect_equal(res$std.dev,apply(TestBootMatrix ,2,sd))
        
})




