context("Test: ecScore() ")


data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

tai.vals <- TAI(PhyloExpressionSetExample)
tdi.vals <- TDI(DivergenceExpressionSetExample)

test_that("ecScore computes correct early conservation scores...",{
        expect_equal(ecScore(tai.vals,
                                           early = 1:2, mid = 3:5, late = 6:7),min(c(mean(tai.vals[3:5] - mean(tai.vals[1:2])), mean(tai.vals[6:7] - mean(tai.vals[1:2])))))
        
        expect_equal(ecScore(tdi.vals,
                             early = 1:2, mid = 3:5, late = 6:7),min(c(mean(tdi.vals[3:5] - mean(tdi.vals[1:2])), mean(tdi.vals[6:7] - mean(tdi.vals[1:2])))))
})
