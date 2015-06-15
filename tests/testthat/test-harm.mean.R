
context("Test: harm.mean() ")

# test function for harm.mean()
test_harmonic_mean <- function(x)
{
        
        N <- length(x)
        ratio <- 1/x
        harm <- 1 / ((1/N) * sum(ratio))
        return(harm)
        
}

test_that("harm.mean() computes correct values",{
        
        expect_equal(harm.mean(c(1,2,4)),12/7)
        expect_equal(harm.mean(1:10),test_harmonic_mean(1:10))
        expect_equal(harm.mean(1:100),test_harmonic_mean(1:100))
        expect_equal(harm.mean(1:1000),test_harmonic_mean(1:1000))
        expect_equal(harm.mean(1:10000),test_harmonic_mean(1:10000))
})

test_that("harm.mean() throws an error when only non-numeric value as passed to the function",{
        
        expect_error(harm.mean("A"),"Please enter a numeric vector.")
        expect_error(harm.mean(as.complex(5)),"Please enter a numeric vector.")
})
