context("Test: PairwiseTest() ")

data(PhyloExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[, 2:9]

test_that(
  "is.ExpressionSet() throws error when no ExpressionSet is entered to PairwiseTest()",
  {
    expect_error(
      PairwiseTest(
        nonStandardExpressionSet,
        modules = list(
          contrast1 = 1:4,
          contrast2 = 5:7
        ),
        altHypothesis = "greater",
        permutations = 1000
      )
    )
  }
)


test_that("p.value is computed..", {
  expect_true(
    PairwiseTest(
      PhyloExpressionSetExample,
      modules = list(
        contrast1 = 1:4,
        contrast2 = 5:7
      ),
      altHypothesis = "greater",
      permutations = 1000
    )$p.value > 0.9
  )
})


test_that("std.dev is computed..", {
  expect_true(length(
    PairwiseTest(
      PhyloExpressionSetExample,
      modules = list(
        contrast1 = 1:4,
        contrast2 = 5:7
      ),
      altHypothesis = "greater",
      permutations = 1000
    )$std.dev
  ) == 7)
})


test_that("lillie.test is NA..", {
  expect_true(is.na(
    PairwiseTest(
      PhyloExpressionSetExample,
      modules = list(
        contrast1 = 1:4,
        contrast2 = 5:7
      ),
      altHypothesis = "greater",
      permutations = 1000
    )$lillie.test
  ))
})

test_that("lillie.test is computed...", {
  
  skip_on_cran()   
  expect_output(PairwiseTest(
    PhyloExpressionSetExample,
    modules = list(
      contrast1 = 1:2,
      contrast2 = 3:7
    ),
    altHypothesis = "greater",
    permutations = 1000,
    lillie.test = TRUE
  )$lillie.test)
  
  
})

test_that("error occurs when module selection does not match number of developmental stages..",
          {
            expect_error(
              PairwiseTest(
                PhyloExpressionSetExample,
                modules = list(
                  contrast1 = 1:2,
                  contrast2 = 3:8
                ),
                altHypothesis = "greater",
                permutations = 1000
              ),
              "The module selection is outside the range of the given ExpressionSet."
              #"The number of stages classified into the two modules does not match the total number of stages stored in the given ExpressionSet."
            )
          })


test_that("error occurs when modules aren't specified...", {
  expect_error(PairwiseTest(PhyloExpressionSetExample,
                                    permutations = 1000))
})



test_that("PairwiseTest() computes correct std.dev and p.values values...",
          {
            skip_on_cran()
            TestBootMatrix <- bootMatrix(PhyloExpressionSetExample, 1000)
            
            res <- PairwiseTest(
              PhyloExpressionSetExample,
              modules = list(
                contrast1 = 1:2,
                contrast2 = 3:7
              ),
              altHypothesis = "greater",
              custom.perm.matrix = TestBootMatrix
            )
            
            estimates <-
              fitdistrplus::fitdist(apply(
                TestBootMatrix,
                1,
                pairScore,
                contrast1 = 1:2,
                contrast2 = 3:7,
                altHypothesis = "greater"
              ),
              distr = "norm")
            
            real_score <- pairScore(TAI(PhyloExpressionSetExample), 1:2, 3:7, altHypothesis = "greater")
            expect_equal(
              res$p.value,
              pnorm(
                real_score,
                mean = estimates$estimate[1],
                sd = estimates$estimate[2],
                lower.tail = FALSE
              )
            )
            expect_equal(res$std.dev, apply(TestBootMatrix , 2, sd))
            
          })
