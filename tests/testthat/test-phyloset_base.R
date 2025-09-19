test_that("BulkPhyloExpressionSet throws error on NA in strata_values", {
    # Create proper data structure: gene_id, strata, expression1, expression2, expression3
    data <- data.frame(
        strata = c(1, 2),
        gene_id = c("Gene1", "Gene2"),
        expr1 = c(1, 4),
        expr2 = c(2, 5),
        expr3 = c(3, 6)
    )
    groups <- c("A", "B", "A")  # 3 groups for 3 expression columns
    legend <- data.frame(stratum = 1:2, label = c("S1", "S2"))
    obj <- as_PhyloExpressionSet(data, groups = groups, name = "Test", species = "Testus", index_type = "TAI", strata_legend = legend)
    expect_error({
        obj@strata_values <- c(NA, 1)
    }, "cannot contain NA values")
})
test_that("PhyloExpressionSet loads correctly", {
    expect_s7_class(example_phyex_set, myTAI::BulkPhyloExpressionSet)
    expect_true(length(example_phyex_set@gene_ids) > 0)
    expect_true(length(example_phyex_set@strata) > 0)
    expect_true(ncol(example_phyex_set@expression) > 0)
    expect_true(nrow(example_phyex_set@expression) > 0)
})

test_that("PhyloExpressionSet properties work", {
    # Test basic properties
    expect_equal(length(example_phyex_set@gene_ids), nrow(example_phyex_set@expression))
    expect_equal(length(example_phyex_set@strata), nrow(example_phyex_set@expression))
    expect_equal(example_phyex_set@num_genes, nrow(example_phyex_set@expression))
    
    # Test computed properties
    expect_true(is.matrix(example_phyex_set@expression_collapsed))
    expect_true(is.data.frame(as_data_frame(example_phyex_set)))
    expect_true(is.numeric(example_phyex_set@TXI))
    expect_true(is.matrix(pTXI(example_phyex_set)))
    
    # Test dimensions
    expect_equal(nrow(pTXI(example_phyex_set)), example_phyex_set@num_genes)
    expect_equal(ncol(pTXI(example_phyex_set)), example_phyex_set@num_identities)
    expect_equal(length(example_phyex_set@TXI), example_phyex_set@num_identities)
})

test_that("PhyloExpressionSet transformation works", {
    # Test log transformation
    log_set <- transform_counts(example_phyex_set, log1p, "log1p")
    expect_s7_class(log_set, myTAI::BulkPhyloExpressionSet)
    expect_equal(nrow(log_set@expression), nrow(example_phyex_set@expression))
    expect_equal(ncol(log_set@expression), ncol(example_phyex_set@expression))
    expect_true(all(log_set@expression >= 0))  # log1p should be non-negative
    
    # Test sqrt transformation
    sqrt_set <- tf(example_phyex_set, sqrt)
    expect_s7_class(sqrt_set, myTAI::BulkPhyloExpressionSet)
    expect_true(all(sqrt_set@expression >= 0))
})

test_that("Gene selection works", {
    # Select subset of genes
    selected_genes <- example_phyex_set@gene_ids[1:10]
    subset_set <- select_genes(example_phyex_set, selected_genes)
    
    expect_s7_class(subset_set, myTAI::BulkPhyloExpressionSet)
    expect_equal(length(subset_set@gene_ids), 10)
    expect_equal(nrow(subset_set@expression), 10)
    expect_true(all(subset_set@gene_ids %in% selected_genes))
})

test_that("Gene removal works", {
    # Remove some genes
    genes_to_remove <- example_phyex_set@gene_ids[1:5]
    filtered_set <- remove_genes(example_phyex_set, genes_to_remove)
    
    expect_s7_class(filtered_set, myTAI::BulkPhyloExpressionSet)
    expect_equal(length(filtered_set@gene_ids), 
                 length(example_phyex_set@gene_ids) - 5)
    expect_true(all(!genes_to_remove %in% filtered_set@gene_ids))
})

test_that("Collapse function works", {
    collapsed_set <- collapse(example_phyex_set)
    expect_s7_class(collapsed_set, myTAI::BulkPhyloExpressionSet)
    expect_equal(nrow(collapsed_set@expression), nrow(example_phyex_set@expression))
    expect_equal(ncol(collapsed_set@expression), example_phyex_set@num_identities)
})

test_that("pTXI calculation works", {
    ptxi_matrix <- pTXI(example_phyex_set)
    expect_true(is.matrix(ptxi_matrix))
    expect_equal(nrow(ptxi_matrix), nrow(example_phyex_set@expression_collapsed))
    expect_equal(ncol(ptxi_matrix), ncol(example_phyex_set@expression_collapsed))
    expect_true(all(is.finite(ptxi_matrix)))
})

test_that("sTXI calculation works", {
    stxi_identity <- sTXI(example_phyex_set, option = "identity")
    stxi_add <- sTXI(example_phyex_set, option = "add")
    
    expect_true(is.matrix(stxi_identity))
    expect_true(is.matrix(stxi_add))
    expect_equal(ncol(stxi_identity), example_phyex_set@num_identities)
    expect_equal(ncol(stxi_add), example_phyex_set@num_identities)
    expect_equal(nrow(stxi_identity), example_phyex_set@num_strata)
    expect_equal(nrow(stxi_add), example_phyex_set@num_strata)
})
