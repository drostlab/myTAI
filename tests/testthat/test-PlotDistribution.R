context("Test: PlotDistribution() ")

data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

nonStandardExpressionSet <- PhyloExpressionSetExample[, 2:9]

test_that(
        "is.ExpressionSet() throughs error when no ExpressionSet is entered to PlotDistribution()",
        {
                expect_error(
                        PlotDistribution(nonStandardExpressionSet)
                        )
        }
)