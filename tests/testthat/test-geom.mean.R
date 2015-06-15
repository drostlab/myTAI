context("Test: geom.mean() ")

# test function for geom.mean()
test_geom_mean <- function(x){
        
        N <- length(x)
        exp(sum(log(x))/N)
        
} 

test_that("geom.mean() computes correct values",{
        
        expect_equal(geom.mean(c(1,4,1/32)),0.5)
        expect_equal(geom.mean(1:10),test_geom_mean(1:10))
        expect_equal(geom.mean(1:100),test_geom_mean(1:100))
        expect_equal(geom.mean(1:1000),test_geom_mean(1:1000))
        expect_equal(geom.mean(1:10000),test_geom_mean(1:10000))
})

test_that("geom.mean() throws an error when only non-numeric value as passed to the function",{
        
        expect_error(geom.mean("A"),"Please enter a numeric vector.")
        expect_error(geom.mean(as.complex(5)),"Please enter a numeric vector.")
})



