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






