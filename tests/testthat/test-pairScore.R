context("Test: pairScore() ")


data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

tai.vals <- TAI(PhyloExpressionSetExample)
tdi.vals <- TDI(DivergenceExpressionSetExample)

test_that("pairScore computes correct pairwise difference scores...",{
  expect_equal(pairScore(tai.vals,
                         contrast1 = 1:2,contrast2 = 3:7,altHypothesis = "greater"),
               mean(tai.vals[1:2]) - mean(tai.vals[3:7]))
  
  expect_equal(pairScore(tdi.vals,
                         contrast1 = 1:2,contrast2 = 3:7, altHypothesis = "greater"),
               mean(tdi.vals[1:2]) - mean(tdi.vals[3:7]))
})

test_that("pairScore computes correct pairwise difference scores...",{
  expect_equal(pairScore(tai.vals,
                         contrast1 = 1:2,contrast2 = 3:7,altHypothesis = "greater"),
               -1*pairScore(tai.vals,
                         contrast1 = 1:2,contrast2 = 3:7,altHypothesis = "less"))
  
  expect_equal(pairScore(tdi.vals,
                         contrast1 = 1:2,contrast2 = 3:7,altHypothesis = "greater"),
               -1*pairScore(tdi.vals,
                            contrast1 = 1:2,contrast2 = 3:7,altHypothesis = "less"))
})