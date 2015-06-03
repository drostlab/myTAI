context("Test: PlotGeneSet() ")

# http://stats.stackexchange.com/questions/118033/best-series-of-colors-to-use-for-differentiating-series-in-publication-quality

equal_df <- function(df1,df2){
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1,df2))
}

data(PhyloExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[ , 2:9] 
example.geneset <- PhyloExpressionSetExample[sample(1:22000,10) , 2]

missing_genes_example <- c("AT4G22950.1","AT4G24540.1", "AT1G69120.1", "AT5G15840.1", "AT4G35900.1", "AT5G10140.1", "AT4G00650.1", "AT1G65480.1", "AT1G22770.1", "AT5G61850.1", "AT2G45660", "AT5G03840")

all_missing_genes_example <- c("AT1G22770", "AT5G61850", "AT2G45660", "AT5G03840")

genes_example <- c("AT4G22950.1","AT4G24540.1", "AT1G69120.1", "AT5G15840.1", "AT4G35900.1", "AT5G10140.1", "AT4G00650.1", "AT1G65480.1", "AT1G22770.1", "AT5G61850.1", "AT2G45660.1", "AT5G03840.1")

TestSubSet <- PhyloExpressionSetExample[match(tolower(genes_example),tolower(PhyloExpressionSetExample[ , 2])), ]
        
test_that("is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotGeneSet()",{
        expect_error(PlotGeneSet(nonStandardExpressionSet,example.geneset),"The present input object does not fulfill the ExpressionSet standard.")
})


test_that("Number of colors and number of input genes are equal when colors are specified..",
          {
                  
                  expect_error(PlotGeneSet(PhyloExpressionSetExample,genes_example, colors = rainbow(10),get.subset = TRUE),"The number of colors and the number of genes do not match.")
                  
          })

test_that("Missing genes are reported by warning....",{
        
        expect_warning(PlotGeneSet(PhyloExpressionSetExample,missing_genes_example,get.subset = TRUE),"Only 10 out of your 12 gene ids could be found in the ExpressionSet.")
        
})

test_that("Error occurs when none of the input gene ids match with the ExpressionSet object",
          {
                  
                  expect_error(PlotGeneSet(PhyloExpressionSetExample,all_missing_genes_example) , "None of your input gene ids could be found in the ExpressionSet.")
                  
          })


test_that("PlotGeneSet() returns the correct subset of genes..",{
        
        expect_true(equal_df(PlotGeneSet(PhyloExpressionSetExample,genes_example,get.subset = TRUE), TestSubSet))
})


