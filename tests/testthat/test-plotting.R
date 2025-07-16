test_that("plot_signature works", {
    # Basic signature plot
    p <- plot_signature(example_phyex_set)
    expect_s3_class(p, "ggplot")
    
    # Test with confidence intervals
    p_ci <- plot_signature(example_phyex_set, show_ci = TRUE)
    expect_s3_class(p_ci, "ggplot")
    
    # Test with different index type
    p_tai <- plot_signature(example_phyex_set, index_type = "TAI")
    expect_s3_class(p_tai, "ggplot")
})

test_that("plot_gene_space works", {
    # Basic gene space plot
    p <- plot_gene_space(example_phyex_set)
    expect_s3_class(p, "ggplot")
    
    # Test with different coloring
    p_strata <- plot_gene_space(example_phyex_set, colour_by = "strata")
    expect_s3_class(p_strata, "ggplot")
    
    # Test with specific genes
    selected_genes <- example_phyex_set@gene_ids[1:20]
    p_genes <- plot_gene_space(example_phyex_set, genes = selected_genes)
    expect_s3_class(p_genes, "ggplot")
})

test_that("plot_gene_profiles works", {
    # Basic gene profiles plot
    selected_genes <- example_phyex_set@gene_ids[1:5]
    p <- plot_gene_profiles(example_phyex_set, genes = selected_genes)
    expect_s3_class(p, "ggplot")
    
    # Test with different transformations
    p_log <- plot_gene_profiles(example_phyex_set, genes = selected_genes, 
                               transformation = "log")
    expect_s3_class(p_log, "ggplot")
    
    p_std_log <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                   transformation = "std_log")
    expect_s3_class(p_std_log, "ggplot")
    
    p_none <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                transformation = "none")
    expect_s3_class(p_none, "ggplot")
    
    # Test with different coloring schemes
    p_strata <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                  colour_by = "strata")
    expect_s3_class(p_strata, "ggplot")
    
    p_stage <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                 colour_by = "stage")
    expect_s3_class(p_stage, "ggplot")
    
    p_manual <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                  colour_by = "manual")
    expect_s3_class(p_manual, "ggplot")
    
    # Test with custom colors
    custom_colors <- c("red", "blue", "green", "purple", "orange")
    names(custom_colors) <- selected_genes
    p_custom <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                  colour_by = "manual", colours = custom_colors)
    expect_s3_class(p_custom, "ggplot")
    
    # Test with show_set_mean
    p_mean <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                show_set_mean = TRUE)
    expect_s3_class(p_mean, "ggplot")
    
    # Test with show_reps
    p_reps <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                show_reps = TRUE)
    expect_s3_class(p_reps, "ggplot")
    
    # Test with faceting by strata
    p_facet <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                 facet_by_strata = TRUE)
    expect_s3_class(p_facet, "ggplot")
    
    # Test with combined options
    p_combined <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                    transformation = "std_log",
                                    colour_by = "strata",
                                    show_set_mean = TRUE,
                                    show_reps = TRUE,
                                    facet_by_strata = TRUE)
    expect_s3_class(p_combined, "ggplot")
    
    # Test without labels
    p_no_labels <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                     show_labels = FALSE)
    expect_s3_class(p_no_labels, "ggplot")
    
    # Test without legend
    p_no_legend <- plot_gene_profiles(example_phyex_set, genes = selected_genes,
                                     show_legend = FALSE)
    expect_s3_class(p_no_legend, "ggplot")
    
    # Test with max_genes parameter (auto gene selection)
    p_auto <- plot_gene_profiles(example_phyex_set, genes = NULL, max_genes = 10)
    expect_s3_class(p_auto, "ggplot")
    
    # Test with larger max_genes
    p_auto_large <- plot_gene_profiles(example_phyex_set, genes = NULL, max_genes = 50)
    expect_s3_class(p_auto_large, "ggplot")
})

test_that("plot_sample_space works", {
    # Test PCA plot
    p_pca <- plot_sample_space(example_phyex_set, method = "PCA")
    expect_s3_class(p_pca, "ggplot")
    
    # Test UMAP plot (if uwot is available)
    if (requireNamespace("uwot", quietly = TRUE)) {
        p_umap <- plot_sample_space(example_phyex_set, method = "UMAP")
        expect_s3_class(p_umap, "ggplot")
    }
})

test_that("plot_mean_var works", {
    p <- plot_mean_var(example_phyex_set)
    expect_s3_class(p, "ggplot")
})

test_that("plot_strata_expression works", {
    # Basic strata expression plot
    p <- plot_strata_expression(example_phyex_set)
    expect_s3_class(p, "ggplot")
    
    # Test with rank version
    p_rank <- plot_strata_expression_rank(example_phyex_set)
    expect_s3_class(p_rank, "ggplot")
})

test_that("plot_distribution_strata works", {
    # Test with a subset of genes
    selected_genes <- example_phyex_set@gene_ids[1:100]
    p <- plot_distribution_strata(example_phyex_set@strata, selected_genes)
    expect_s3_class(p, "ggplot")
})

test_that("plot_contribution works", {
    p <- plot_contribution(example_phyex_set)
    expect_s3_class(p, "ggplot")
})

test_that("plot_distribution_expression works", {
    # Basic distribution plot
    p <- plot_distribution_expression(example_phyex_set)
    expect_s3_class(p, "ggplot")
    
    # Test with different parameters
    p_strata <- plot_distribution_expression(example_phyex_set, show_strata = TRUE)
    expect_s3_class(p_strata, "ggplot")
    
    p_no_conditions <- plot_distribution_expression(example_phyex_set, show_conditions = FALSE)
    expect_s3_class(p_no_conditions, "ggplot")
})

test_that("Plotting functions handle edge cases", {
    # Test with single gene
    single_gene <- example_phyex_set@gene_ids[1]
    p_single <- plot_gene_profiles(example_phyex_set, genes = single_gene)
    expect_s3_class(p_single, "ggplot")
    
    # Test with single gene and all options
    p_single_full <- plot_gene_profiles(example_phyex_set, genes = single_gene,
                                       transformation = "std_log",
                                       colour_by = "manual",
                                       show_set_mean = TRUE,
                                       show_reps = TRUE,
                                       show_labels = TRUE,
                                       show_legend = TRUE,
                                       facet_by_strata = TRUE)
    expect_s3_class(p_single_full, "ggplot")
    
    # Test with empty gene list (should use auto selection)
    p_empty <- plot_gene_profiles(example_phyex_set, genes = NULL, max_genes = 5)
    expect_s3_class(p_empty, "ggplot")
    
    # Test with very large max_genes (should be capped)
    p_large <- plot_gene_profiles(example_phyex_set, genes = NULL, max_genes = 10000)
    expect_s3_class(p_large, "ggplot")
    
    # Test with max_genes = 1
    p_max1 <- plot_gene_profiles(example_phyex_set, genes = NULL, max_genes = 1)
    expect_s3_class(p_max1, "ggplot")
    
    # Test different combinations of show parameters
    p_combo1 <- plot_gene_profiles(example_phyex_set, genes = example_phyex_set@gene_ids[1:3],
                                  show_set_mean = TRUE, show_reps = FALSE, show_labels = FALSE)
    expect_s3_class(p_combo1, "ggplot")
    
    p_combo2 <- plot_gene_profiles(example_phyex_set, genes = example_phyex_set@gene_ids[1:3],
                                  show_set_mean = FALSE, show_reps = TRUE, show_labels = TRUE)
    expect_s3_class(p_combo2, "ggplot")
    
    # Test with colour_by = "stage" and faceting
    p_stage_facet <- plot_gene_profiles(example_phyex_set, genes = example_phyex_set@gene_ids[1:5],
                                       colour_by = "stage", facet_by_strata = TRUE)
    expect_s3_class(p_stage_facet, "ggplot")
    
    # Test with colour_by = "manual" and faceting
    p_manual_facet <- plot_gene_profiles(example_phyex_set, genes = example_phyex_set@gene_ids[1:5],
                                        colour_by = "manual", facet_by_strata = TRUE)
    expect_s3_class(p_manual_facet, "ggplot")
    
    # Test with different subset sizes
    p_subset2 <- plot_gene_profiles(example_phyex_set, genes = example_phyex_set@gene_ids[1:2])
    expect_s3_class(p_subset2, "ggplot")
    
    p_subset10 <- plot_gene_profiles(example_phyex_set, genes = example_phyex_set@gene_ids[1:10])
    expect_s3_class(p_subset10, "ggplot")
    
    # Test with very small gene set for signature plot
    small_set <- select_genes(example_phyex_set, example_phyex_set@gene_ids[1:5])
    p_small <- plot_signature(small_set)
    expect_s3_class(p_small, "ggplot")
})

test_that("plot_signature_transformed works", {
    # Test with different transformations
    transforms <- list(
        identity = function(x) x,
        log1p = function(x) log1p(x)
    )
    
    p <- plot_signature_transformed(example_phyex_set, transforms)
    expect_s3_class(p, "ggplot")
})

test_that("plot_signature_multiple works", {
    # Create multiple datasets for comparison
    log_set <- transform_counts(example_phyex_set, log1p, "log1p")
    sqrt_set <- transform_counts(example_phyex_set, sqrt, "sqrt")
    
    datasets <- list(
        "Original" = example_phyex_set,
        "Log" = log_set,
        "Sqrt" = sqrt_set
    )
    
    p <- plot_signature_multiple(datasets)
    expect_s3_class(p, "ggplot")
})
