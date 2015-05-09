context("Test: CollapseReplicates() ")

# adapted from: https://github.com/hadley/dplyr/blob/master/tests/testthat/test-arrange.r 
equal_df <- function(df1, df2) {
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1, df2))
}

data(PhyloExpressionSetExample)

# test df: mean, nrep = 2, stage.names = c("S1","S2","S3")
TestCollapsedDF1 <- cbind(PhyloExpressionSetExample[1:5, 1:2],apply(PhyloExpressionSetExample[1:5, 3:4],1,mean),apply(PhyloExpressionSetExample[1:5, 5:6],1,mean),apply(PhyloExpressionSetExample[1:5, 7:8],1,mean))
colnames(TestCollapsedDF1)[3:5] <- c("S1","S2","S3")

# test df: mean, nrep = c(2,2,3), stage.names = c("S1","S2","S3")
TestCollapsedDF2 <- cbind(PhyloExpressionSetExample[1:5, 1:2],apply(PhyloExpressionSetExample[1:5, 3:4],1,mean),apply(PhyloExpressionSetExample[1:5, 5:6],1,mean),apply(PhyloExpressionSetExample[1:5, 7:9],1,mean))
colnames(TestCollapsedDF2)[3:5] <- c("S1","S2","S3")

# test df: median, nrep = 2, stage.names = c("S1","S2","S3")
TestCollapsedDF3 <- cbind(PhyloExpressionSetExample[1:5, 1:2],apply(PhyloExpressionSetExample[1:5, 3:4],1,median),apply(PhyloExpressionSetExample[1:5, 5:6],1,median),apply(PhyloExpressionSetExample[1:5, 7:8],1,median))
colnames(TestCollapsedDF3)[3:5] <- c("S1","S2","S3")

# test df: max, nrep = 2, stage.names = c("S1","S2","S3")
TestCollapsedDF4 <- cbind(PhyloExpressionSetExample[1:5, 1:2],apply(PhyloExpressionSetExample[1:5, 3:4],1,max),apply(PhyloExpressionSetExample[1:5, 5:6],1,max),apply(PhyloExpressionSetExample[1:5, 7:8],1,max))
colnames(TestCollapsedDF4)[3:5] <- c("S1","S2","S3")


# test df: median, nrep = c(2,2,3), stage.names = c("S1","S2","S3")
TestCollapsedDF5 <- cbind(PhyloExpressionSetExample[1:5, 1:2],apply(PhyloExpressionSetExample[1:5, 3:4],1,median),apply(PhyloExpressionSetExample[1:5, 5:6],1,median),apply(PhyloExpressionSetExample[1:5, 7:9],1,median))
colnames(TestCollapsedDF5)[3:5] <- c("S1","S2","S3")

# test df: max, nrep = c(2,2,3), stage.names = c("S1","S2","S3")
TestCollapsedDF6 <- cbind(PhyloExpressionSetExample[1:5, 1:2],apply(PhyloExpressionSetExample[1:5, 3:4],1,max),apply(PhyloExpressionSetExample[1:5, 5:6],1,max),apply(PhyloExpressionSetExample[1:5, 7:9],1,max))
colnames(TestCollapsedDF6)[3:5] <- c("S1","S2","S3")


test_that("CollapseReplicates() works with const. nrep = 2 and FUN = mean", {
        
        expect_true(equal_df(CollapseReplicates(PhyloExpressionSetExample[1:5,1:8],2,mean,c("S1","S2","S3")), TestCollapsedDF1))
        
})

test_that("CollapseReplicates() works with const. nrep = 2 and FUN = 'mean' (quoted)", {
        
        expect_true(equal_df(CollapseReplicates(PhyloExpressionSetExample[1:5,1:8],2,"mean",c("S1","S2","S3")), TestCollapsedDF1))
        
})


test_that("CollapseReplicates() works with const. nrep = 2 and FUN = median", {
        
        expect_true(equal_df(CollapseReplicates(PhyloExpressionSetExample[1:5,1:8],2,median,c("S1","S2","S3")), TestCollapsedDF3))
        
})

test_that("CollapseReplicates() works with const. nrep = 2 and FUN = max", {
        
        expect_true(equal_df(CollapseReplicates(PhyloExpressionSetExample[1:5,1:8],2,max,c("S1","S2","S3")), TestCollapsedDF4))
        
})

test_that("CollapseReplicates() works with variable nrep = c(2,2,3) and FUN = mean", {
        
        expect_true(equal_df(CollapseReplicates(PhyloExpressionSetExample[1:5,1:9],c(2,2,3),mean,c("S1","S2","S3")), TestCollapsedDF2))
        
})


test_that("CollapseReplicates() works with variable nrep = c(2,2,3) and FUN = median", {
        
        expect_true(equal_df(CollapseReplicates(PhyloExpressionSetExample[1:5,1:9],c(2,2,3),median,c("S1","S2","S3")), TestCollapsedDF5))
        
})

test_that("CollapseReplicates() works with variable nrep = c(2,2,3) and FUN = max", {
        
        expect_true(equal_df(CollapseReplicates(PhyloExpressionSetExample[1:5,1:9],c(2,2,3),max,c("S1","S2","S3")), TestCollapsedDF6))
        
})

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to CollapseReplicates()",{
        expect_error(CollapseReplicates(nonStandardExpressionSet,2,mean),"The present input object does not fulfill the ExpressionSet standard.")
})


TestCollapsedDF_stage.names <- TestCollapsedDF1
colnames(TestCollapsedDF_stage.names)[3:5] <- paste0("X",1:3)
        
test_that("default stage.names work properly when stage.names = NULL",{
        
        expect_true(equal_df(CollapseReplicates(PhyloExpressionSetExample[1:5,1:8],2,mean), TestCollapsedDF_stage.names))
        
})

