context("Test : tfStability() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

equal_df <- function(df1, df2) {
  rownames(df1) <- NULL
  rownames(df2) <- NULL
  isTRUE(all.equal(df1, df2))
}

nonStandardExpressionSet <- PhyloExpressionSetExample[, 2:9]

test_that("is.ExpressionSet() throws error when no ExpressionSet is entered to tfStability()",
          {
            expect_error(
              tfStability(nonStandardExpressionSet)
            )
            
          })

test_that("tfStability() throws error when an unavailable TestStatistic is entered to tfStability()",
          {
            expect_error(
              tfStability(PhyloExpressionSetExample, TestStatistic = "WavyLineTest")
            )
            
          })

test_that("tfStability() throws error when modules are not specified",
          {
            expect_error(
              tfStability(PhyloExpressionSetExample, TestStatistic = "ReductiveHourglassTest")
            )
            
          })

test_that("error occurs when module selection does not match number of developmental stages..",
          {
            expect_error(
              tfStability(PhyloExpressionSetExample, 
                          TestStatistic = "ReductiveHourglassTest",
                          modules = list(
                            early = 1:2,
                            mid = 3:5,
                            late = 6:8))
            )
            
          })

test_that("error occurs when transformation is not available",
          {
            expect_error(
              tfStability(PhyloExpressionSetExample, 
                          TestStatistic = "ReductiveHourglassTest",
                          transforms = c("rank", "sqwuared"),
                          modules = list(
                            early = 1:2,
                            mid = 3:5,
                            late = 6:7),
                          permutations = 100)
            )
            
          })

test_that("tfStability() computes correct p-values ", {
  expect_true(tfStability(ExpressionSet = PhyloExpressionSetExample,
                TestStatistic = "ReductiveHourglassTest",
                transforms = c("log2", "sqrt", "none"),
                modules = list(early = 1:2, mid = 3:5, late = 6:7),
                permutations = 100)[2] < 0.05
  )
})