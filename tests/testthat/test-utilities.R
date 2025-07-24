test_that("PS_colours works", {
    # Test basic color generation
    colors_5 <- PS_colours(5)
    expect_true(is.character(colors_5))
    expect_equal(length(colors_5), 5)
    expect_true(all(grepl("^#", colors_5)))  # Should be hex colors
    
    # Test with different numbers
    colors_10 <- PS_colours(10)
    expect_equal(length(colors_10), 10)
    
    colors_1 <- PS_colours(1)
    expect_equal(length(colors_1), 1)
})

test_that("omit_matrix works", {
    # Test omit matrix calculation
    omit_mat <- omit_matrix(example_phyex_set)
    expect_true(is.matrix(omit_mat))
    expect_equal(nrow(omit_mat), length(example_phyex_set@gene_ids))
    expect_equal(ncol(omit_mat), example_phyex_set@num_identities)
    expect_true(all(is.finite(omit_mat)))
    
    # Check row names
    expect_true(all(grepl("^\\(-\\)", rownames(omit_mat))))
})

test_that("age.apply works", {
    # Test with column means
    result_means <- age.apply(example_phyex_set, colMeans)
    expect_true(is.matrix(result_means))
    expect_equal(nrow(result_means), example_phyex_set@num_strata)
    expect_equal(ncol(result_means), example_phyex_set@num_identities)
    
    # Test with as.list = TRUE
    result_list <- age.apply(example_phyex_set, colMeans, as.list = TRUE)
    expect_true(is.list(result_list))
    expect_equal(length(result_list), example_phyex_set@num_strata)
    
    # Test with custom function
    result_var <- age.apply(example_phyex_set, function(x) apply(x, 2, var))
    expect_true(is.matrix(result_var))
})

test_that("genes_filter_dynamic works", {
    # Test basic filtering
    expr_matrix <- log1p(example_phyex_set@expression_collapsed)
    filtered <- genes_filter_dynamic(expr_matrix, thr = 0.8)
    expect_true(is.matrix(filtered))
    expect_true(nrow(filtered) <= nrow(expr_matrix))
    expect_equal(ncol(filtered), ncol(expr_matrix))
    
    # Test with different thresholds
    filtered_90 <- genes_filter_dynamic(expr_matrix, thr = 0.9)
    filtered_50 <- genes_filter_dynamic(expr_matrix, thr = 0.5)
    expect_true(nrow(filtered_90) <= nrow(filtered_50))
})

test_that(".to_std_expr works", {
    # Test standardization (internal function)
    expr_matrix <- log1p(example_phyex_set@expression_collapsed)
    std_expr <- myTAI:::.to_std_expr(expr_matrix)
    expect_true(is.matrix(std_expr))
    expect_equal(dim(std_expr), dim(expr_matrix))
    
    # Check that genes with non-zero variance are standardized
    valid_genes <- apply(expr_matrix, 1, function(x) stats::var(x) > 0)
    if (any(valid_genes)) {
        gene_means <- rowMeans(std_expr[valid_genes, , drop = FALSE])
        gene_sds <- apply(std_expr[valid_genes, , drop = FALSE], 1, stats::sd)
        expect_true(all(abs(gene_means) < 1e-10))  # Should be ~0
        expect_true(all(abs(gene_sds - 1) < 1e-10))  # Should be ~1
    }
})

test_that("get_angles works", {
    # Test angle calculation (internal function)
    expr_matrix <- log1p(example_phyex_set@expression_collapsed)
    std_expr <- myTAI:::.to_std_expr(expr_matrix)
    angles <- myTAI:::get_angles(std_expr)
    
    expect_true(is.numeric(angles))
    expect_equal(length(angles), nrow(std_expr))
    expect_true(all(angles >= -pi & angles <= pi))
})

test_that("Gene pattern functions work", {
    S <- 6  # Number of stages
    
    # Test pattern functions
    early_pattern <- myTAI:::early_gene(S)
    mid_pattern <- myTAI:::mid_gene(S)
    late_pattern <- myTAI:::late_gene(S)
    rev_mid_pattern <- myTAI:::rev_mid_gene(S)
    
    expect_equal(length(early_pattern), S)
    expect_equal(length(mid_pattern), S)
    expect_equal(length(late_pattern), S)
    expect_equal(length(rev_mid_pattern), S)
    
    # Test relationships
    expect_equal(late_pattern, -early_pattern)
    expect_equal(rev_mid_pattern, -mid_pattern)
    
    # Test mod_pi
    angles <- c(-4, -pi, 0, pi, 4)
    wrapped <- mod_pi(angles)
    expect_true(all(wrapped >= -pi & wrapped <= pi))
})

test_that("tfPS works", {
    # Test phylostratum transformation
    tf_set <- tfPS(example_phyex_set, transform = "qr")
    expect_s7_class(tf_set, PhyloExpressionSet)
    expect_equal(length(tf_set@strata), length(example_phyex_set@strata))
    # Check that the strata are factors with the same level names
    expect_true(is.factor(tf_set@strata))
    expect_true(all(levels(tf_set@strata) %in% levels(example_phyex_set@strata)))
    
    # Test with quantilerank
    tf_set_qr <- tfPS(example_phyex_set, transform = "quantilerank")
    expect_s7_class(tf_set_qr, PhyloExpressionSet)
    expect_equal(tf_set@strata, tf_set_qr@strata)
})

test_that("Data conversion functions work", {
    # Test match_map
    # Create dummy phylomap
    phylomap <- data.frame(
        Stratum = example_phyex_set@strata[1:100],
        GeneID = example_phyex_set@gene_ids[1:100]
    )
    
    # Create dummy expression data
    expr_data <- data.frame(
        GeneID = example_phyex_set@gene_ids[1:100],
        example_phyex_set@expression_collapsed[1:100, ]
    )
    
    matched_set <- match_map(expr_data, phylomap)
    expect_s7_class(matched_set, myTAI::BulkPhyloExpressionSet)
    expect_equal(nrow(matched_set@expression), 100)
})
