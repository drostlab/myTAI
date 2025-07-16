test_that("PhyloExpressionSet loads correctly", {
    expect_s7_class(example_phyex_set, myTAI::PhyloExpressionSet)
    expect_true(length(example_phyex_set@gene_ids) > 0)
    expect_true(length(example_phyex_set@strata) > 0)
    expect_true(ncol(example_phyex_set@counts) > 0)
    expect_true(nrow(example_phyex_set@counts) > 0)
})

test_that("PhyloExpressionSet properties work", {
    # Test basic properties
    expect_equal(length(example_phyex_set@gene_ids), nrow(example_phyex_set@counts))
    expect_equal(length(example_phyex_set@strata), nrow(example_phyex_set@counts))
    expect_equal(example_phyex_set@num_genes, nrow(example_phyex_set@counts))
    
    # Test computed properties
    expect_true(is.matrix(example_phyex_set@counts_collapsed))
    expect_true(is.data.frame(example_phyex_set@data))
    expect_true(is.numeric(example_phyex_set@TXI))
    expect_true(is.matrix(example_phyex_set@pTXI))
    
    # Test dimensions
    expect_equal(nrow(example_phyex_set@pTXI), example_phyex_set@num_genes)
    expect_equal(ncol(example_phyex_set@pTXI), example_phyex_set@num_conditions)
    expect_equal(length(example_phyex_set@TXI), example_phyex_set@num_conditions)
})

test_that("PhyloExpressionSet transformation works", {
    # Test log transformation
    log_set <- transform_counts(example_phyex_set, log1p, "log1p")
    expect_s7_class(log_set, myTAI::PhyloExpressionSet)
    expect_equal(nrow(log_set@counts), nrow(example_phyex_set@counts))
    expect_equal(ncol(log_set@counts), ncol(example_phyex_set@counts))
    expect_true(all(log_set@counts >= 0))  # log1p should be non-negative
    
    # Test sqrt transformation
    sqrt_set <- tf(example_phyex_set, sqrt)
    expect_s7_class(sqrt_set, myTAI::PhyloExpressionSet)
    expect_true(all(sqrt_set@counts >= 0))
})

test_that("Gene selection works", {
    # Select subset of genes
    selected_genes <- example_phyex_set@gene_ids[1:10]
    subset_set <- select_genes(example_phyex_set, selected_genes)
    
    expect_s7_class(subset_set, myTAI::PhyloExpressionSet)
    expect_equal(length(subset_set@gene_ids), 10)
    expect_equal(nrow(subset_set@counts), 10)
    expect_true(all(subset_set@gene_ids %in% selected_genes))
})

test_that("Gene removal works", {
    # Remove some genes
    genes_to_remove <- example_phyex_set@gene_ids[1:5]
    filtered_set <- remove_genes(example_phyex_set, genes_to_remove)
    
    expect_s7_class(filtered_set, myTAI::PhyloExpressionSet)
    expect_equal(length(filtered_set@gene_ids), 
                 length(example_phyex_set@gene_ids) - 5)
    expect_true(all(!genes_to_remove %in% filtered_set@gene_ids))
})

test_that("Collapse function works", {
    collapsed_set <- collapse(example_phyex_set)
    expect_s7_class(collapsed_set, myTAI::PhyloExpressionSet)
    expect_equal(nrow(collapsed_set@counts), nrow(example_phyex_set@counts))
    expect_equal(ncol(collapsed_set@counts), example_phyex_set@num_conditions)
})

test_that("pTXI calculation works", {
    ptxi_matrix <- pTXI(example_phyex_set@counts_collapsed, example_phyex_set@strata)
    expect_true(is.matrix(ptxi_matrix))
    expect_equal(nrow(ptxi_matrix), nrow(example_phyex_set@counts_collapsed))
    expect_equal(ncol(ptxi_matrix), ncol(example_phyex_set@counts_collapsed))
    expect_true(all(is.finite(ptxi_matrix)))
})

test_that("sTXI calculation works", {
    stxi_identity <- sTXI(example_phyex_set, option = "identity")
    stxi_add <- sTXI(example_phyex_set, option = "add")
    
    expect_true(is.matrix(stxi_identity))
    expect_true(is.matrix(stxi_add))
    expect_equal(ncol(stxi_identity), example_phyex_set@num_conditions)
    expect_equal(ncol(stxi_add), example_phyex_set@num_conditions)
    expect_equal(nrow(stxi_identity), example_phyex_set@num_strata)
    expect_equal(nrow(stxi_add), example_phyex_set@num_strata)
})
