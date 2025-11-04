test_that("genes_top_variance works", {
    # Test default parameters
    top_var_genes <- genes_top_variance(example_phyex_set)
    expect_true(is.character(top_var_genes))
    expect_true(length(top_var_genes) > 0)
    expect_true(all(top_var_genes %in% example_phyex_set@gene_ids))
    
    # Test with different percentile
    top_var_90 <- genes_top_variance(example_phyex_set, top_p = 0.9)
    top_var_95 <- genes_top_variance(example_phyex_set, top_p = 0.95)
    expect_true(length(top_var_90) >= length(top_var_95))
    
    # Test with very low percentile
    top_var_50 <- genes_top_variance(example_phyex_set, top_p = 0.5)
    expect_true(length(top_var_50) >= length(top_var_90))
})

test_that("genes_top_mean works", {
    # Test default parameters
    top_expr_genes <- genes_top_mean(example_phyex_set)
    expect_true(is.character(top_expr_genes))
    expect_true(length(top_expr_genes) > 0)
    expect_true(all(top_expr_genes %in% example_phyex_set@gene_ids))
    
    # Test with different percentile
    top_expr_90 <- genes_top_mean(example_phyex_set, top_p = 0.9)
    top_expr_95 <- genes_top_mean(example_phyex_set, top_p = 0.95)
    expect_true(length(top_expr_90) >= length(top_expr_95))
})

test_that("genes_lowly_expressed works", {
    # Test default threshold
    lowly_expr_genes <- genes_lowly_expressed(example_phyex_set)
    expect_true(is.character(lowly_expr_genes))
    expect_true(all(lowly_expr_genes %in% example_phyex_set@gene_ids))
    
    # Test with different threshold
    lowly_expr_5 <- genes_lowly_expressed(example_phyex_set, threshold = 5)
    lowly_expr_10 <- genes_lowly_expressed(example_phyex_set, threshold = 10)
    expect_true(length(lowly_expr_5) <= length(lowly_expr_10))
    
    # Test with very high threshold (should return most genes)
    lowly_expr_high <- genes_lowly_expressed(example_phyex_set, threshold = 1000)
    expect_true(length(lowly_expr_high) >= length(lowly_expr_5))
})

test_that("Gene selection functions return valid genes", {
    # Test that selected genes actually exist in the dataset
    top_var <- genes_top_variance(example_phyex_set, top_p = 0.8)
    top_expr <- genes_top_mean(example_phyex_set, top_p = 0.8)
    lowly_expr <- genes_lowly_expressed(example_phyex_set, threshold = 5)
    
    # All should be valid gene IDs
    expect_true(all(top_var %in% example_phyex_set@gene_ids))
    expect_true(all(top_expr %in% example_phyex_set@gene_ids))
    expect_true(all(lowly_expr %in% example_phyex_set@gene_ids))
    
    # No duplicates
    expect_equal(length(top_var), length(unique(top_var)))
    expect_equal(length(top_expr), length(unique(top_expr)))
    expect_equal(length(lowly_expr), length(unique(lowly_expr)))
})

test_that("Gene selection with extreme parameters", {
    # Test with top_p = 0 (should return all genes >= 0th percentile)
    top_var_0 <- genes_top_variance(example_phyex_set, top_p = 0)
    expect_equal(length(top_var_0), length(example_phyex_set@gene_ids))  # Should be all genes
    
    # Test with top_p = 0.99 (should return top 1% of genes)
    top_var_99 <- genes_top_variance(example_phyex_set, top_p = 0.99)
    expect_true(length(top_var_99) <= 1000)  # Should be small subset
    
    # Test with threshold = 0 (should return very few genes)
    lowly_0 <- genes_lowly_expressed(example_phyex_set, threshold = 0)
    expect_true(length(lowly_0) <= length(example_phyex_set@gene_ids))
})

test_that("Gene selection with top_k parameter works", {
    # Test with top_k for variance
    top_var_100 <- genes_top_variance(example_phyex_set, top_k = 100)
    expect_equal(length(top_var_100), 100)
    expect_true(all(top_var_100 %in% example_phyex_set@gene_ids))
    
    # Test with top_k for mean
    top_mean_50 <- genes_top_mean(example_phyex_set, top_k = 50)
    expect_equal(length(top_mean_50), 50)
    expect_true(all(top_mean_50 %in% example_phyex_set@gene_ids))
    
    # Test that top_k takes precedence over top_p
    top_var_k <- genes_top_variance(example_phyex_set, top_p = 0.5, top_k = 10)
    expect_equal(length(top_var_k), 10)
    
    # Test with top_k = 0
    top_var_0k <- genes_top_variance(example_phyex_set, top_k = 0)
    expect_equal(length(top_var_0k), 0)
    
    # Test with top_k larger than total genes
    n_genes <- length(example_phyex_set@gene_ids)
    top_var_all <- genes_top_variance(example_phyex_set, top_k = n_genes + 100)
    expect_equal(length(top_var_all), n_genes)
})
