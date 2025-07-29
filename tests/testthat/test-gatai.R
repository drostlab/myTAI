test_that("destroy_pattern works with gataiR", {
    # Skip if gataiR is not available
    skip_if_not_installed("gataiR")
    
    # Run GATAI with minimal parameters for testing
    gatai_result <- destroy_pattern(example_phyex_set_old, 
                                    num_runs = 3, 
                                    max_generations = 100,
                                    extended_analysis = TRUE,
                                    min_pval = 0.0,
                                    always_return_genes = TRUE)
    
    # Check that the result has the expected structure
    expect_type(gatai_result, "list")
    expect_type(gatai_result$removed_genes, "character")
    expect_type(gatai_result$runs, "list")
    expect_length(gatai_result$runs, 3)  # Should have 3 runs
    
    # Check that some genes were removed
    expect_gt(length(gatai_result$removed_genes), 0)

    # Create a temporary directory for testing
    temp_dir <- tempdir()
    test_dir <- file.path(temp_dir, "gatai_test")
    
    # Save PDF
    pdf_path <- save_gatai_results_pdf(example_phyex_set_old, gatai_result, test_dir,
                                       prefix = "test_gatai")
    
    # Check that the PDF was created
    expect_true(file.exists(pdf_path))
    expect_match(pdf_path, "test_gatai\\.pdf$")
    
    # Check that the file has some content (basic check)
    file_info <- file.info(pdf_path)
    expect_gt(file_info$size, 1000)  # Should be at least 1KB
    
    # Clean up
    unlink(test_dir, recursive = TRUE)

})

test_that("destroy_pattern handles empty output from gataiR", {
    # Skip if gataiR is not available
    skip_if_not_installed("gataiR")
    res <- destroy_pattern(example_phyex_set_old, analysis_dir="gatai", max_generations = 5, num_runs = 3)
    expect_equal(length(res$removed_genes), 0)
})