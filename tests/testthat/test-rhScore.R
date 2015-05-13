context("Test: rhScore() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)


data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

tai.vals <- TAI(PhyloExpressionSetExample)
tdi.vals <- TDI(DivergenceExpressionSetExample)

test_that("rhScore computes correct reductive hourglass scores...",{
        expect_equal(rhScore(tai.vals,
                             early = 1:2, mid = 3:5, late = 6:7,method = "min", scoringMethod = "mean-mean"),min(c(mean(tai.vals[1:2] - mean(tai.vals[3:5])), mean(tai.vals[6:7] - mean(tai.vals[3:5])))))
        
        expect_equal(rhScore(tdi.vals,
                             early = 1:2, mid = 3:5, late = 6:7,method = "min", scoringMethod = "mean-mean"),min(c(mean(tdi.vals[1:2] - mean(tdi.vals[3:5])), mean(tdi.vals[6:7] - mean(tdi.vals[3:5])))))
})


