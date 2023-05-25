context("Test : tf() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

equal_df <- function(df1, df2) {
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1, df2))
}

nonStandardExpressionSet <- PhyloExpressionSetExample[, 2:9]

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to tf()",
          {
                  expect_error(
                          tf(nonStandardExpressionSet)
                  )
                  
          })

Test.tf.Func <- function(ExpressionSet, FUN) {
  f <- match.fun(FUN)
  ncols <- ncol(ExpressionSet)
  res <-
    cbind(ExpressionSet[, 1:2], apply(ExpressionSet[, 3:ncols], 2, f))
  return(res)
}


test_that("tf() computes correct transformation values ", {
        expect_true(equal_df(
                tf(PhyloExpressionSetExample , log2),
                Test.tf.Func(PhyloExpressionSetExample, log2)
        ))
        expect_true(equal_df(
                tf(DivergenceExpressionSetExample , log2),
                Test.tf.Func(DivergenceExpressionSetExample, log2)
        ))
        
        expect_true(equal_df(
                tf(PhyloExpressionSetExample , sqrt),
                Test.tf.Func(PhyloExpressionSetExample, sqrt)
        ))
        expect_true(equal_df(
                tf(DivergenceExpressionSetExample , sqrt),
                Test.tf.Func(DivergenceExpressionSetExample, sqrt)
        ))
        
        expect_true(equal_df(
                tf(PhyloExpressionSetExample , log2),
                Test.tf.Func(PhyloExpressionSetExample, log2)
        ))
        expect_true(equal_df(
                tf(DivergenceExpressionSetExample , sqrt),
                Test.tf.Func(DivergenceExpressionSetExample, sqrt)
        ))
        
        expect_true(equal_df(
                tf(PhyloExpressionSetExample , function(x)
                        x / 2),
                Test.tf.Func(PhyloExpressionSetExample, function(x)
                        x / 2)
        ))
        expect_true(equal_df(
                tf(DivergenceExpressionSetExample , function(x)
                        x / 2),
                Test.tf.Func(DivergenceExpressionSetExample, function(x)
                        x / 2)
        ))
})

# For the newer functionalities

Test.tf.Func <- function(ExpressionSet, FUN, pseudocount = 0, integerise = FALSE) {
  f <- match.fun(FUN)
  
  ExpressionMatrix <- as.matrix(ExpressionSet[ , -c(1,2)] + pseudocount)

  if(integerise){
    ExpressionMatrix <- round(ExpressionMatrix, digits = 0)
  }
  
  f <- match.fun(FUN)
  res_mat <- f(ExpressionMatrix)
  
  res <- base::cbind(ExpressionSet[ , c(1,2)], base::as.data.frame(res_mat))
  return(res)
}

test_that("tf() computes correct transformation values ", {
  expect_true(equal_df(
    tf(PhyloExpressionSetExample , log2, 
       pseudocount = 1, integerise = TRUE),
    Test.tf.Func(PhyloExpressionSetExample, log2, 
                 pseudocount = 1, integerise = TRUE)
  ))
  expect_true(equal_df(
    tf(DivergenceExpressionSetExample , log2, 
       pseudocount = 1, integerise = TRUE),
    Test.tf.Func(DivergenceExpressionSetExample, log2, 
                 pseudocount = 1, integerise = TRUE)
  ))
  
  expect_true(equal_df(
    tf(PhyloExpressionSetExample , 
       FUN = function(x) apply(x, 2, base::rank), 
       pseudocount = 1, integerise = TRUE),
    Test.tf.Func(PhyloExpressionSetExample, 
                 FUN = function(x) apply(x, 2, base::rank), 
                 pseudocount = 1, integerise = TRUE)
  ))
  expect_true(equal_df(
    tf(DivergenceExpressionSetExample , 
       FUN = function(x) apply(x, 2, base::rank), 
       pseudocount = 1, integerise = TRUE),
    Test.tf.Func(DivergenceExpressionSetExample, 
                 FUN = function(x) apply(x, 2, base::rank), 
                 pseudocount = 1, integerise = TRUE)
  ))
  })

test_that("tf() throws error when pseudocount is not a numeric value",
          {
            expect_error(
              tf(PhyloExpressionSetExample, pseudocount = "A")
            )
            
          })

test_that("tf() throws error when pseudocount is not a single value",
          {
            expect_error(
              tf(PhyloExpressionSetExample, pseudocount = c(1,2))
            )
            
          })