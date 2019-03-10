context("Test: PlotEnrichment() ")

# http://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists
# http://en.wikipedia.org/wiki/Fisher%27s_exact_test



equal_df <- function(df1, df2) {
        rownames(df1) <- NULL
        rownames(df2) <- NULL
        isTRUE(all.equal(df1, df2))
}

missing_genes_example <-
        c(
                "AT4G22950.1",
                "AT4G24540.1",
                "AT1G69120.1",
                "AT5G15840.1",
                "AT4G35900.1",
                "AT5G10140.1",
                "AT4G00650.1",
                "AT1G65480.1",
                "AT1G22770.1",
                "AT5G61850.1",
                "AT2G45660",
                "AT5G03840"
        )

all_missing_genes_example <-
        c("AT1G22770", "AT5G61850", "AT2G45660", "AT5G03840")

genes_example <-
        c(
                "AT4G22950.1",
                "AT4G24540.1",
                "AT1G69120.1",
                "AT5G15840.1",
                "AT4G35900.1",
                "AT5G10140.1",
                "AT4G00650.1",
                "AT1G65480.1",
                "AT1G22770.1",
                "AT5G61850.1",
                "AT2G45660.1",
                "AT5G03840.1"
        )

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

set.seed(123)
test_set <- sample(PhyloExpressionSetExample[, 2], 22000)


nonStandardExpressionSet <- PhyloExpressionSetExample[, 2:9]

test_that(
        "is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotEnrichment()",
        {
                expect_error(
                        PlotEnrichment(nonStandardExpressionSet, test_set, legendName = "PS")
                        )
        }
)

test_that(
        "error occurs when test.set inlcudes more genes than are available in the
        input ExpressionSet",
        {
                expect_error(
                        PlotEnrichment(PhyloExpressionSetExample[sample(1:22000, 19000), ], test_set, legendName = "PS") ,
                        "Your input GeneID vector stores more elements than are available in your ExpressionSet object..."
                )
                
        }
)


test_that("Missing genes are reported by warning....", {
        expect_warning(
                PlotEnrichment(
                        PhyloExpressionSetExample,
                        missing_genes_example,
                        legendName = "PS",
                        plot.bars = FALSE
                ),
                "Only 10 out of your 12 gene ids could be found in the ExpressionSet."
        )
        
})

test_that("Error occurs when none of the input gene ids match with the ExpressionSet object",
          {
                  expect_error(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  all_missing_genes_example,
                                  legendName = "PS"
                          ) ,
                          "None of your input gene ids could be found in the ExpressionSet."
                  )
                  
          })



TestSetPS1_3 <-
        PhyloExpressionSetExample[which(PhyloExpressionSetExample[, 1] %in% c(1, 2, 3)), 2]

TestSetPS1_2 <-
        PhyloExpressionSetExample[which(PhyloExpressionSetExample[, 1] %in% c(1, 2)), 2]

TestSetPS1 <-
        PhyloExpressionSetExample[which(PhyloExpressionSetExample[, 1] %in% c(1)), 2]


test_that("PlotEnrichment() executes with only 3 out of 12 PS present in the test.set",
          {
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1_3,
                                  measure = "log-foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[4:12 , 2] == 0
                  )))
                  
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1_3,
                                  measure = "log-foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[1:3 , 2] > 0
                  )))
                  
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1_3,
                                  measure = "foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[4:12 , 2] == 0
                  )))
                  
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1_3,
                                  measure = "foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[1:3 , 2] > 0
                  )))
          })


test_that("PlotEnrichment() executes with only 2 out of 12 PS present in the test.set",
          {
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1_2,
                                  measure = "log-foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[4:12 , 2] == 0
                  )))
                  
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1_2,
                                  measure = "log-foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[1:3 , 2] > 0
                  )))
                  
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1_2,
                                  measure = "foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[4:12 , 2] == 0
                  )))
                  
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1_2,
                                  measure = "foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[1:3 , 2] > 0
                  )))
          })


test_that("PlotEnrichment() executes with only 1 out of 12 PS present in the test.set",
          {
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1,
                                  measure = "log-foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[4:12 , 2] == 0
                  )))
                  
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1,
                                  measure = "log-foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[1:3 , 2] > 0
                  )))
                  
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1,
                                  measure = "foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[4:12 , 2] == 0
                  )))
                  
                  expect_true(all(which(
                          PlotEnrichment(
                                  PhyloExpressionSetExample,
                                  TestSetPS1,
                                  measure = "foldchange",
                                  legendName = "PS",
                                  plot.bars = FALSE
                          )$enrichment.matrix[1:3 , 2] > 0
                  )))
          })


# construct an example PhyloExpressionSet

PES.Test <- data.frame(
        PS = c(
                rep(1, 4),
                rep(2, 2),
                rep(3, 6),
                rep(4, 1),
                rep(5, 12),
                rep(6, 3),
                rep(7, 3)
        ),
        GeneID = paste0("AT", 1:31),
        stringsAsFactors = FALSE
)


