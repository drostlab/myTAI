# Test single-cell functionality with S7 PhyloExpressionSet objects

test_that("ScPhyloExpressionSet can be created", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    expect_s7_class(example_phyex_set_sc, myTAI::ScPhyloExpressionSet)
    expect_true(length(example_phyex_set_sc@gene_ids) > 0)
    expect_true(length(example_phyex_set_sc@strata) > 0)
    expect_s4_class(example_phyex_set_sc@seurat, "Seurat")
    expect_true(is.factor(example_phyex_set_sc@groups))
})

test_that("ScPhyloExpressionSet properties work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test basic properties
    expect_equal(length(example_phyex_set_sc@gene_ids), example_phyex_set_sc@num_genes)
    expect_equal(length(example_phyex_set_sc@strata), example_phyex_set_sc@num_genes)
    expect_true(is.numeric(example_phyex_set_sc@TXI_sample))
    expect_true(is.matrix(example_phyex_set_sc@expression_collapsed))
})

test_that("Single-cell plotting functions work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test signature plot for single-cell
    p_sig <- plot_signature(example_phyex_set_sc)
    expect_s3_class(p_sig, "ggplot")
    
    # Test with individual cells
    p_cells <- plot_signature(example_phyex_set_sc, show_reps = TRUE)
    expect_s3_class(p_cells, "ggplot")
    
    # Test gene heatmap
    p_heat <- plot_gene_heatmap(example_phyex_set_sc, top_p = 0.1)
    expect_s3_class(p_heat, "ggplot")
    
    # Test with cell sampling
    p_heat_cells <- plot_gene_heatmap(example_phyex_set_sc, reps = TRUE, max_cells_per_type = 3)
    expect_s3_class(p_heat_cells, "ggplot")
})




test_that("Single-cell transformations work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test transformation
    log_sc_set <- transform_counts(example_phyex_set_sc, log1p, "log1p")
    expect_s7_class(log_sc_set, myTAI::ScPhyloExpressionSet)
    expect_true(grepl("log1p", log_sc_set@name))
    
    # Test gene selection
    selected_genes <- example_phyex_set_sc@gene_ids[1:10]
    subset_sc_set <- select_genes(example_phyex_set_sc, selected_genes)
    expect_equal(length(subset_sc_set@gene_ids), 10)
    expect_s7_class(subset_sc_set, myTAI::ScPhyloExpressionSet)
})

test_that("Single-cell sample space plots work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test PCA plot (should use existing Seurat reduction)
    p_pca <- plot_sample_space(example_phyex_set_sc, method = "PCA")
    expect_s3_class(p_pca, "ggplot")
    
    # Test UMAP plot (should use existing Seurat reduction)
    p_umap <- plot_sample_space(example_phyex_set_sc, method = "UMAP")
    expect_s3_class(p_umap, "ggplot")
    
    # Test with different coloring
    p_color <- plot_sample_space(example_phyex_set_sc, method = "PCA", colour_by = "identity")
    expect_s3_class(p_color, "ggplot")
})

test_that("Single-cell gene space plots work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Basic gene space plot
    p <- plot_gene_space(example_phyex_set_sc)
    expect_s3_class(p, "ggplot")
    
    # Test with specific genes
    selected_genes <- example_phyex_set_sc@gene_ids[1:20]
    p_genes <- plot_gene_space(example_phyex_set_sc, genes = selected_genes)
    expect_s3_class(p_genes, "ggplot")
})

test_that("Single-cell gene profiles work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Basic gene profiles plot for single-cell
    selected_genes <- example_phyex_set_sc@gene_ids[1:5]
    p <- plot_gene_profiles(example_phyex_set_sc, genes = selected_genes)
    expect_s3_class(p, "ggplot")
    
    # Test with replicates (individual cells)
    p_reps <- plot_gene_profiles(example_phyex_set_sc, genes = selected_genes, show_reps = TRUE)
    expect_s3_class(p_reps, "ggplot")
})

test_that("Single-cell distribution plots work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test expression distribution plots
    p_dist <- plot_distribution_expression(example_phyex_set_sc)
    expect_s3_class(p_dist, "ggplot")
    
    # Test with strata distributions
    p_strata <- plot_distribution_expression(example_phyex_set_sc, show_strata = TRUE)
    expect_s3_class(p_strata, "ggplot")
})

test_that("Single-cell specific parameters work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test max_cells_per_type parameter in heatmap
    p_heat1 <- plot_gene_heatmap(example_phyex_set_sc, reps = TRUE, max_cells_per_type = 2)
    expect_s3_class(p_heat1, "ggplot")
    
    p_heat2 <- plot_gene_heatmap(example_phyex_set_sc, reps = TRUE, max_cells_per_type = 5)
    expect_s3_class(p_heat2, "ggplot")
    
    # Test that single-cell signature plots use violin plots
    p_sig <- plot_signature(example_phyex_set_sc)
    expect_s3_class(p_sig, "ggplot")
    
    # Test single-cell signature with individual cells
    p_sig_cells <- plot_signature(example_phyex_set_sc, show_reps = TRUE)
    expect_s3_class(p_sig_cells, "ggplot")
})

test_that("Single-cell mean-var and strata plots work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test mean-variance plot
    p_mv <- plot_mean_var(example_phyex_set_sc)
    expect_s3_class(p_mv, "ggplot")
    
    # Test with gene highlighting
    selected_genes <- example_phyex_set_sc@gene_ids[1:5]
    p_mv_highlight <- plot_mean_var(example_phyex_set_sc, highlight_genes = selected_genes)
    expect_s3_class(p_mv_highlight, "ggplot")
    
    # Test strata expression plot
    p_strata <- plot_strata_expression(example_phyex_set_sc)
    expect_s3_class(p_strata, "ggplot")
})

test_that("Single-cell contribution plots work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test phylostratum contribution plot
    p_contrib <- plot_contribution(example_phyex_set_sc)
    expect_s3_class(p_contrib, "ggplot")
})

test_that("Single-cell multiple signature plots work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Create transformed datasets for comparison
    log_sc_set <- transform_counts(example_phyex_set_sc, log1p, "log1p")
    sqrt_sc_set <- transform_counts(example_phyex_set_sc, sqrt, "sqrt")
    
    datasets <- list(
        "Original" = example_phyex_set_sc,
        "Log" = log_sc_set,
        "Sqrt" = sqrt_sc_set
    )
    
    p_multi <- plot_signature_multiple(datasets)
    expect_s3_class(p_multi, "ggplot")
})

test_that("ScPhyloExpressionSet properties are computed correctly", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test that TXI_sample has correct length (number of cells)
    expect_equal(length(example_phyex_set_sc@TXI_sample), 
                 length(example_phyex_set_sc@groups))
    
    # Test that expression_collapsed has correct dimensions
    expect_equal(nrow(example_phyex_set_sc@expression_collapsed), 
                 example_phyex_set_sc@num_genes)
    expect_equal(ncol(example_phyex_set_sc@expression_collapsed), 
                 example_phyex_set_sc@num_identities)
    
    # Test that TXI has correct length (number of cell types)
    expect_equal(length(example_phyex_set_sc@TXI), 
                 example_phyex_set_sc@num_identities)
})

test_that("ScPhyloExpressionSet edge cases work", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    

    
    # Test with small gene subset
    small_genes <- example_phyex_set_sc@gene_ids[1:5]
    small_sc_set <- select_genes(example_phyex_set_sc, small_genes)
    
    p_small_heat <- plot_gene_heatmap(small_sc_set, genes = small_genes)
    expect_s3_class(p_small_heat, "ggplot")
})

test_that("plot_signature with flexible identities works", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test basic plot with default identity
    p1 <- plot_signature(example_phyex_set_sc)
    expect_s3_class(p1, "ggplot")
    
    # Test with primary identity specified
    p2 <- plot_signature(example_phyex_set_sc, primary_identity = "day")
    expect_s3_class(p2, "ggplot")
    
    # Test with both primary and secondary identity (coloring)
    p3 <- plot_signature(example_phyex_set_sc, 
                        primary_identity = "day", 
                        secondary_identity = "condition")
    expect_s3_class(p3, "ggplot")
    
    # Test with faceting by secondary identity
    p4 <- plot_signature(example_phyex_set_sc, 
                        primary_identity = "day", 
                        secondary_identity = "condition",
                        facet_by_secondary = TRUE)
    expect_s3_class(p4, "ggplot")
    
    # Test with show_reps = FALSE
    p5 <- plot_signature(example_phyex_set_sc, 
                        primary_identity = "batch",
                        show_reps = FALSE)
    expect_s3_class(p5, "ggplot")
})

test_that("plot_signature identity validation works", {
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggforce")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test error for invalid primary identity
    expect_error(
        plot_signature(example_phyex_set_sc, primary_identity = "nonexistent"),
        "Primary identity 'nonexistent' not found in metadata"
    )
    
    # Test error for invalid secondary identity
    expect_error(
        plot_signature(example_phyex_set_sc, 
                      primary_identity = "day",
                      secondary_identity = "nonexistent"),
        "Secondary identity 'nonexistent' not found in metadata"
    )
})

test_that("available_identities function works", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    identities <- available_identities(example_phyex_set_sc)
    expect_type(identities, "character")
    expect_true(length(identities) > 0)
    expect_true("groups" %in% identities)
    expect_true("day" %in% identities)
    expect_true("condition" %in% identities)
    expect_true("batch" %in% identities)
})

test_that("set_identity_colours function works", {
    skip_if_not_installed("Seurat")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test setting colours for groups identity
    colours <- c("TypeA" = "red", "TypeB" = "blue", "TypeC" = "green")
    sc_set_coloured <- set_identity_colours(example_phyex_set_sc, "groups", colours)
    
    expect_s7_class(sc_set_coloured, myTAI::ScPhyloExpressionSet)
    expect_true(length(sc_set_coloured@identity_colours) > 0)
    expect_equal(sc_set_coloured@identity_colours[["groups"]], colours)
    
    # Test setting colours for day identity
    day_colours <- c("Day1" = "lightblue", "Day3" = "blue", "Day5" = "darkblue", "Day7" = "navy")
    sc_set_day_coloured <- set_identity_colours(example_phyex_set_sc, "day", day_colours)
    
    expect_s7_class(sc_set_day_coloured, myTAI::ScPhyloExpressionSet)
    expect_equal(sc_set_day_coloured@identity_colours[["day"]], day_colours)
    
    # Test error for invalid identity
    expect_error(
        set_identity_colours(example_phyex_set_sc, "nonexistent", colours),
        "Identity 'nonexistent' not found in metadata"
    )
    
    # Test error for invalid colours format
    expect_error(
        set_identity_colours(example_phyex_set_sc, "groups", c("red", "blue", "green")),
        "Colours must be a named character vector"
    )
})

test_that("plot_signature shows color message when using defaults", {
    skip_if_not_installed("Seurat")
    skip_if_not_installed("ggforce")
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test that message is shown when using default colors
    expect_message(
        plot_signature(example_phyex_set_sc, primary_identity = "day"),
        "Using default colors for identity 'day'"
    )
    
    # Test that no message is shown when using custom color parameter
    expect_silent(
        plot_signature(example_phyex_set_sc, primary_identity = "day", colour = "red")
    )
    
    # Test that no message is shown when using custom colors from object
    colors <- c("Day1" = "red", "Day3" = "blue", "Day5" = "green", "Day7" = "purple")
    sc_set_colored <- set_identity_colours(example_phyex_set_sc, "day", colors)
    expect_silent(
        plot_signature(sc_set_colored, primary_identity = "day")
    )
})
