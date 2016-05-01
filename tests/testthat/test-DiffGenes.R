context("Test: DiffGenes() ")

data(PhyloExpressionSetExample)

equal_df <- function(df1, df2) {
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1, df2))
}

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 

test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to DiffGenes()",{
        expect_error(DiffGenes(ExpressionSet = nonStandardExpressionSet,
                               nrep          = 2,
                               comparison    = "below",
                               method        = "foldchange",
                               stage.names   = c("S1","S2","S3"))
                     , "The present input object does not fulfill the ExpressionSet standard.")
})



test_that("Correct fold change values are computed..",{
        
        set.seed(123)
        ExampleMicroarrayTable <- PhyloExpressionSetExample[sample(1:22000,10) , 1:8] 
        
        stage.1 <- apply( cbind(ExampleMicroarrayTable[ , 3],ExampleMicroarrayTable[ , 4]) , 1 , myTAI::geom.mean )
        stage.2 <- apply( cbind(ExampleMicroarrayTable[ , 5],ExampleMicroarrayTable[ , 6]) , 1 , myTAI::geom.mean )
        
        
        expect_equal(DiffGenes(ExpressionSet = ExampleMicroarrayTable,
                               nrep          = 2,
                               comparison    = "below",
                               method        = "foldchange",
                               stage.names   = c("S1","S2","S3"))[ , 3], stage.1 / stage.2)
        
})




test_that("Correct log fold change values are computed..",{
        
        set.seed(123)
        ExampleMicroarrayTable <- PhyloExpressionSetExample[sample(1:22000,10) , 1:8] 
        
        stage.1 <- apply( cbind(ExampleMicroarrayTable[ , 3],ExampleMicroarrayTable[ , 4]) , 1 , myTAI::geom.mean )
        stage.2 <- apply( cbind(ExampleMicroarrayTable[ , 5],ExampleMicroarrayTable[ , 6]) , 1 , myTAI::geom.mean )
        
        
        expect_equal(DiffGenes(ExpressionSet = tf(ExampleMicroarrayTable,log2),
                               nrep          = 2,
                               comparison    = "below",
                               method        = "log-foldchange",
                               stage.names   = c("S1","S2","S3"))[ , 3], log2(stage.1) - log2(stage.2))
        
})

test_that("Correct p-values based on the t-test are computed..",{
        
        set.seed(123)
        ExampleMicroarrayTable <- PhyloExpressionSetExample[sample(1:22000,10) , 1:8] 
        
        stage.comp <- apply( cbind(ExampleMicroarrayTable[ , 3:4],ExampleMicroarrayTable[ , 5:6]) , 1 , function(x) t.test(x[1:2],x[3:4])$p.value)
        
        
        expect_equal(DiffGenes(ExpressionSet = ExampleMicroarrayTable,
                               nrep          = 2,
                               comparison    = "below",
                               method        = "t.test",
                               stage.names   = c("S1","S2","S3"))[ , 3], as.numeric(stage.comp))
        
})


test_that("Correct p-values based on the wilcox.test are computed..",{

        set.seed(123)
        ExampleMicroarrayTable <- PhyloExpressionSetExample[sample(1:22000,10) , 1:8] 
        
        stage.comp <- apply( cbind(ExampleMicroarrayTable[ , 3:4],ExampleMicroarrayTable[ , 5:6]) , 1 , function(x) wilcox.test(x[1:2],x[3:4])$p.value)
        
        
        expect_equal(DiffGenes(ExpressionSet = ExampleMicroarrayTable,
                               nrep          = 2,
                               comparison    = "below",
                               method        = "wilcox.test",
                               stage.names   = c("S1","S2","S3"))[ , 3], as.numeric(stage.comp))
        
})



