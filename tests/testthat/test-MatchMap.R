context("Test: MatchMap() ")

data(PhyloExpressionSetExample)

# adapted from: https://github.com/hadley/dplyr/blob/master/tests/testthat/test-arrange.r 
equal_df <- function(df1, df2) {
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1, df2))
}

TestMap <- PhyloExpressionSetExample[c(243,2456,6891,14819,21205), 1:2]
TestExpressionSet <- PhyloExpressionSetExample[ , 2:9]

test_that("MatchMap result matches dplyr::inner_join ",{
        MM <- MatchMap(TestMap,TestExpressionSet)
        IJ <- dplyr::inner_join(TestMap,TestExpressionSet,by = "GeneID")
        
        expect_true(equal_df(MM[order(MM[ , "GeneID"]), c(1,3:9)],
                             IJ[order(IJ[ , "GeneID"]), c(1,3:9)]))
        expect_true(isTRUE(all.equal.character(as.character(MM[order(MM[ , "GeneID"]), 2]),as.character(IJ[order(IJ[ , "GeneID"]), 2]))))
        
})
            
test_that("error occurs when duplicate gene ids are stored in either the PhyloMap or ExpressionSet",{
        TestMap2 <- PhyloExpressionSetExample[c(243,2456,243,6891,14819,2456,21205,14819), 1:2]
        
        expect_error(MatchMap(TestMap2,TestExpressionSet),"You have duplicate Gene IDs in your Map. Please enter only unique Gene IDs.")
        
        expect_error(MatchMap(TestMap,PhyloExpressionSetExample[c(243,2456,243,6891,14819,2456,21205,14819), 2:9]),"You have duplicate Gene IDs in your ExpressionMatrix. Please enter only unique Gene IDs, or specify the 'accumulate' argument.")
        
        
})


test_that("MatchMap works properly when duplicates are removed via remove.duplicates = TRUE",{
        
        TestMap3 <- PhyloExpressionSetExample[c(243,2456,243,6891,14819,21205,14819), 1:2]
        
        MM <- MatchMap(TestMap3,TestExpressionSet,remove.duplicates = TRUE)
        IJ <- dplyr::inner_join(TestMap3[-c(1,5), ],TestExpressionSet,by = "GeneID")
        
        expect_true(equal_df(MM[order(MM[ , "GeneID"]), c(1,3:9)],
                             IJ[order(IJ[ , "GeneID"]), c(1,3:9)]))
        expect_true(isTRUE(all.equal.character(as.character(MM[order(MM[ , "GeneID"]), 2]),as.character(IJ[order(IJ[ , "GeneID"]), 2]))))
        
}) 

test_that("MatchMap works properly when duplicates in ExpressionMatrix are accumulated via accumulate = max",{
        
        MM <- MatchMap(TestMap,PhyloExpressionSetExample[c(243,2456,243,6891,14819,21205,14819), 2:9], accumulate = max)
        IJ <- dplyr::inner_join(TestMap,TestExpressionSet,by = "GeneID")
        
        expect_true(equal_df(MM[order(MM[ , "GeneID"]), c(1,3:9)],
                             IJ[order(IJ[ , "GeneID"]), c(1,3:9)]))
        expect_true(isTRUE(all.equal.character(as.character(MM[order(MM[ , "GeneID"]), 2]),as.character(IJ[order(IJ[ , "GeneID"]), 2]))))
        
}) 




            
            
            
            
            
            