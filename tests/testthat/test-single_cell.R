# Test single-cell functionality with S7 PhyloExpressionSet objects

test_that("ScPhyloExpressionSet can be created", {
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    expect_s7_class(example_phyex_set_sc, myTAI::ScPhyloExpressionSet)
    expect_true(length(example_phyex_set_sc@gene_ids) > 0)
    expect_true(length(example_phyex_set_sc@strata) > 0)
    expect_true(is.data.frame(example_phyex_set_sc@metadata))
    expect_true(is.factor(example_phyex_set_sc@groups))
})

test_that("ScPhyloExpressionSet properties work", {
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test basic properties
    expect_equal(length(example_phyex_set_sc@gene_ids), example_phyex_set_sc@num_genes)
    expect_equal(length(example_phyex_set_sc@strata), example_phyex_set_sc@num_genes)
    expect_true(is.numeric(example_phyex_set_sc@TXI_sample))
    expect_true(is.matrix(example_phyex_set_sc@expression_collapsed))
})

test_that("Single-cell plotting functions work", {
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
    p_heat_cells <- plot_gene_heatmap(example_phyex_set_sc, show_reps = TRUE, max_cells_per_type = 3)
    expect_s3_class(p_heat_cells, "ggplot")
})




test_that("Single-cell transformations work", {
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
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test PCA plot (should use stored reductions)
    p_pca <- plot_sample_space(example_phyex_set_sc, method = "PCA")
    expect_s3_class(p_pca, "ggplot")
    
    # Test UMAP plot (should use stored reductions)
    p_umap <- plot_sample_space(example_phyex_set_sc, method = "UMAP")
    expect_s3_class(p_umap, "ggplot")
    
    # Test with different coloring
    p_color <- plot_sample_space(example_phyex_set_sc, method = "PCA", colour_by = "identity")
    expect_s3_class(p_color, "ggplot")
})

test_that("Single-cell gene space plots work", {
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
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test max_cells_per_type parameter in heatmap
    p_heat1 <- plot_gene_heatmap(example_phyex_set_sc, show_reps = TRUE, max_cells_per_type = 2)
    expect_s3_class(p_heat1, "ggplot")
    
    p_heat2 <- plot_gene_heatmap(example_phyex_set_sc, show_reps = TRUE, max_cells_per_type = 5)
    expect_s3_class(p_heat2, "ggplot")
    
    # Test that single-cell signature plots use violin plots
    p_sig <- plot_signature(example_phyex_set_sc)
    expect_s3_class(p_sig, "ggplot")
    
    # Test single-cell signature with individual cells
    p_sig_cells <- plot_signature(example_phyex_set_sc, show_reps = TRUE)
    expect_s3_class(p_sig_cells, "ggplot")
})

test_that("Single-cell mean-var and strata plots work", {
    
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
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test phylostratum contribution plot
    p_contrib <- plot_contribution(example_phyex_set_sc)
    expect_s3_class(p_contrib, "ggplot")
})

test_that("Single-cell multiple signature plots work", {
    
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
    
    # Test that TXI_sample values are reasonable (not all NA, finite values)
    expect_true(any(!is.na(example_phyex_set_sc@TXI_sample)))
    finite_txi <- example_phyex_set_sc@TXI_sample[!is.na(example_phyex_set_sc@TXI_sample)]
    expect_true(all(is.finite(finite_txi)))
    expect_true(all(finite_txi > 0))  # TXI values should be positive
})

test_that("C++ TXI implementation works with ScPhyloExpressionSet", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test that C++ TXI calculation works
    txi_values <- example_phyex_set_sc@TXI_sample
    expect_type(txi_values, "double")
    expect_equal(length(txi_values), example_phyex_set_sc@num_samples)
    
    # Test with different batch sizes if we have enough cells
    if (example_phyex_set_sc@num_samples > 100) {
        # Test direct C++ function call
        txi_cpp_small_batch <- cpp_txi_sc(
            example_phyex_set_sc@expression, 
            example_phyex_set_sc@strata_values, 
            batch_size = 50,
            ncores = 1
        )
        
        txi_cpp_large_batch <- cpp_txi_sc(
            example_phyex_set_sc@expression, 
            example_phyex_set_sc@strata_values, 
            batch_size = 200,
            ncores = 1
        )
        
        # Results should be identical regardless of batch size
        expect_equal(txi_cpp_small_batch, txi_cpp_large_batch, tolerance = 1e-12)
    }
})

test_that("ScPhyloExpressionSet edge cases work", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test with small gene subset
    small_genes <- example_phyex_set_sc@gene_ids[1:5]
    small_sc_set <- select_genes(example_phyex_set_sc, small_genes)
    
    p_small_heat <- plot_gene_heatmap(small_sc_set, genes = small_genes)
    expect_s3_class(p_small_heat, "ggplot")
})

test_that("plot_signature with flexible identities works", {
    
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



test_that("idents_colours property works", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test setting colours for groups identity
    colours <- c("TypeA" = "red", "TypeB" = "blue", "TypeC" = "green")
    example_phyex_set_sc@idents_colours[["groups"]] <- colours
    
    expect_s7_class(example_phyex_set_sc, myTAI::ScPhyloExpressionSet)
    expect_true(length(example_phyex_set_sc@idents_colours) > 0)
    expect_equal(example_phyex_set_sc@idents_colours[["groups"]], colours)
    
    # Test setting colours for day identity
    day_colours <- c("Day1" = "lightblue", "Day3" = "blue", "Day5" = "darkblue", "Day7" = "navy")
    example_phyex_set_sc@idents_colours[["day"]] <- day_colours
    
    expect_s7_class(example_phyex_set_sc, myTAI::ScPhyloExpressionSet)
    expect_equal(example_phyex_set_sc@idents_colours[["day"]], day_colours)
})

test_that("plot_signature shows color message when using defaults", {
    
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
    example_phyex_set_sc@idents_colours[["day"]] <- colors
    expect_silent(
        plot_signature(example_phyex_set_sc, primary_identity = "day")
    )
})

test_that("ScPhyloExpressionSet metadata and groups consistency", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test that groups property matches selected_idents column in metadata
    expected_groups <- example_phyex_set_sc@metadata[[example_phyex_set_sc@selected_idents]]
    expect_equal(as.character(example_phyex_set_sc@groups), as.character(expected_groups))
    
    # Test that metadata has the correct number of rows (one per cell)
    expect_equal(nrow(example_phyex_set_sc@metadata), example_phyex_set_sc@num_samples)
    
    # Test that selected_idents points to a valid column
    expect_true(example_phyex_set_sc@selected_idents %in% colnames(example_phyex_set_sc@metadata))
})

test_that("ScPhyloExpressionSet reductions property works", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test that reductions is a list
    expect_type(example_phyex_set_sc@reductions, "list")
    
    # Test that stored reductions have correct dimensions (cells x components)
    if (length(example_phyex_set_sc@reductions) > 0) {
        for (reduction_name in names(example_phyex_set_sc@reductions)) {
            reduction_data <- example_phyex_set_sc@reductions[[reduction_name]]
            expect_true(is.matrix(reduction_data) || is.data.frame(reduction_data))
            expect_equal(nrow(reduction_data), example_phyex_set_sc@num_samples)
        }
    }
})

test_that("downsample functions work correctly", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test downsample_expression function
    downsampled_expr <- downsample_expression(
        example_phyex_set_sc@expression, 
        example_phyex_set_sc@groups, 
        downsample = 3
    )
    
    expect_true(is.matrix(downsampled_expr))
    expect_equal(nrow(downsampled_expr), example_phyex_set_sc@num_genes)
    expect_true(ncol(downsampled_expr) <= length(example_phyex_set_sc@groups))
    expect_true(all(!is.na(colnames(downsampled_expr))))
    
    # Test downsample method for ScPhyloExpressionSet
    downsampled_set <- downsample(example_phyex_set_sc, downsample = 2)
    expect_s7_class(downsampled_set, myTAI::ScPhyloExpressionSet)
    expect_true(downsampled_set@num_samples <= example_phyex_set_sc@num_samples)
    expect_equal(downsampled_set@num_genes, example_phyex_set_sc@num_genes)
    
    # Test that metadata is properly filtered
    expect_equal(nrow(downsampled_set@metadata), downsampled_set@num_samples)
    
    # Test that reductions are properly filtered if they exist
    if (length(example_phyex_set_sc@reductions) > 0) {
        for (reduction_name in names(example_phyex_set_sc@reductions)) {
            if (reduction_name %in% names(downsampled_set@reductions)) {
                expect_equal(nrow(downsampled_set@reductions[[reduction_name]]), 
                           downsampled_set@num_samples)
            }
        }
    }
})

test_that("ScPhyloExpressionSet cached properties work", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test that pseudobulk expression is cached correctly
    expr_collapsed1 <- example_phyex_set_sc@expression_collapsed
    expr_collapsed2 <- example_phyex_set_sc@expression_collapsed
    
    # Should be identical (cached)
    expect_identical(expr_collapsed1, expr_collapsed2)
    
    # Test that TXI values are cached correctly
    txi1 <- example_phyex_set_sc@TXI
    txi2 <- example_phyex_set_sc@TXI
    expect_identical(txi1, txi2)
    
    # Test that sample-level TXI values are cached correctly
    txi_sample1 <- example_phyex_set_sc@TXI_sample
    txi_sample2 <- example_phyex_set_sc@TXI_sample
    expect_identical(txi_sample1, txi_sample2)
})

test_that("ScPhyloExpressionSet metadata setter validation works", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Create a copy to test on
    test_set <- example_phyex_set_sc
    
    # Test that metadata setter validates row count
    invalid_metadata <- example_phyex_set_sc@metadata[1:5, ]  # Too few rows
    expect_error(
        {test_set@metadata <- invalid_metadata}
    )
    
    # Test that metadata setter converts character columns to factors
    new_metadata <- example_phyex_set_sc@metadata
    new_metadata$new_char_col <- as.character(new_metadata[[1]])  # Add character column
    test_set@metadata <- new_metadata
    
    expect_true(is.factor(test_set@metadata$new_char_col))
})

test_that("ScPhyloExpressionSet selected_idents setter works", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Create a copy to test on
    test_set <- example_phyex_set_sc
    original_selected <- test_set@selected_idents
    
    # Test switching to a different identity
    available_cols <- colnames(test_set@metadata)
    alternative_col <- available_cols[available_cols != original_selected][1]
    
    test_set@selected_idents <- alternative_col
    expect_equal(test_set@selected_idents, alternative_col)
    
    # Test that groups property updates accordingly
    expected_groups <- test_set@metadata[[alternative_col]]
    expect_equal(as.character(test_set@groups), as.character(expected_groups))
})

test_that("Expression matrix consistency checks", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test that expression matrix has correct dimensions
    expect_equal(nrow(example_phyex_set_sc@expression), example_phyex_set_sc@num_genes)
    expect_equal(ncol(example_phyex_set_sc@expression), example_phyex_set_sc@num_samples)
    
    # Test that expression matrix column names match cell identifiers
    if (!is.null(colnames(example_phyex_set_sc@expression))) {
        expect_equal(length(colnames(example_phyex_set_sc@expression)), 
                    example_phyex_set_sc@num_samples)
    }
    
    # Test that expression matrix row names match gene_ids
    if (!is.null(rownames(example_phyex_set_sc@expression))) {
        expect_equal(rownames(example_phyex_set_sc@expression), 
                    example_phyex_set_sc@gene_ids)
    }
    
    # Test that expression matrix is sparse (if expected)
    expect_true(inherits(example_phyex_set_sc@expression, "Matrix") || 
                is.matrix(example_phyex_set_sc@expression))
})

test_that("Strata consistency checks", {
    
    skip_if(is.null(tryCatch(example_phyex_set_sc, error = function(e) NULL)), 
            "example_phyex_set_sc not available")
    
    # Test strata_values property
    expect_equal(length(example_phyex_set_sc@strata_values), 
                example_phyex_set_sc@num_genes)
    expect_true(is.numeric(example_phyex_set_sc@strata_values))
    
    # Test num_strata consistency
    expect_equal(example_phyex_set_sc@num_strata, 
                length(unique(example_phyex_set_sc@strata_values)))
})

test_that("ScPhyloExpressionSet_from_seurat constructor works", {
    skip_if_not_installed("Seurat")
    
    # Create a minimal Seurat object for testing
    set.seed(42)
    n_genes <- 100
    n_cells <- 50
    
    counts <- Matrix::Matrix(
        stats::rnbinom(n_genes * n_cells, size = 3, mu = 50),
        nrow = n_genes,
        ncol = n_cells,
        sparse = TRUE
    )
    rownames(counts) <- paste0("Gene-", 1:n_genes)
    colnames(counts) <- paste0("Cell-", 1:n_cells)
    
    # Create metadata
    metadata <- data.frame(
        cell_type = factor(sample(c("TypeA", "TypeB"), n_cells, replace = TRUE)),
        condition = factor(sample(c("Ctrl", "Treat"), n_cells, replace = TRUE)),
        row.names = colnames(counts)
    )
    
    # Create Seurat object
    seurat_obj <- Seurat::CreateSeuratObject(
        counts = counts,
        meta.data = metadata,
        project = "TestSeurat"
    )
    Seurat::Idents(seurat_obj) <- "cell_type"
    
    # Add some basic processing to get reductions
    suppressWarnings({
        seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
        seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = 20, verbose = FALSE)
        seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
        seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = 5, verbose = FALSE)
    })
    
    # Create strata
    strata <- factor(sample(1:5, n_genes, replace = TRUE))
    names(strata) <- rownames(counts)
    
    # Test constructor
    sc_set <- ScPhyloExpressionSet_from_seurat(
        seurat = seurat_obj,
        strata = strata,
        name = "Test Seurat Constructor"
    )
    
    # Test that object was created correctly
    expect_s7_class(sc_set, myTAI::ScPhyloExpressionSet)
    expect_equal(sc_set@num_genes, n_genes)
    expect_equal(sc_set@num_samples, n_cells)
    expect_equal(sc_set@name, "Test Seurat Constructor")
    
    # Test that metadata was extracted correctly
    expect_true(is.data.frame(sc_set@metadata))
    expect_true("cell_type" %in% colnames(sc_set@metadata))
    expect_true("condition" %in% colnames(sc_set@metadata))
    expect_equal(nrow(sc_set@metadata), n_cells)
    
    # Test that active identities were used correctly
    expect_equal(sc_set@selected_idents, "active_ident")
    
    # Test that reductions were extracted
    expect_type(sc_set@reductions, "list")
    if (length(sc_set@reductions) > 0) {
        expect_true("pca" %in% names(sc_set@reductions) || "PCA" %in% names(sc_set@reductions))
    }
    
    # Test with custom selected_idents
    sc_set_custom <- ScPhyloExpressionSet_from_seurat(
        seurat = seurat_obj,
        strata = strata,
        selected_idents = "condition",
        name = "Test Custom Idents"
    )
    
    expect_equal(sc_set_custom@selected_idents, "condition")
    expected_groups <- sc_set_custom@metadata[["condition"]]
    expect_equal(as.character(sc_set_custom@groups), as.character(expected_groups))
})

test_that("match_map_sc_seurat function works", {
    skip_if_not_installed("Seurat")
    
    # Create a minimal Seurat object for testing
    set.seed(123)
    n_genes <- 50
    n_cells <- 30
    
    counts <- Matrix::Matrix(
        stats::rnbinom(n_genes * n_cells, size = 3, mu = 30),
        nrow = n_genes,
        ncol = n_cells,
        sparse = TRUE
    )
    rownames(counts) <- paste0("TestGene-", 1:n_genes)
    colnames(counts) <- paste0("TestCell-", 1:n_cells)
    # Create Seurat object
    seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = "TestMap")
    
    # Create phylomap (only include some genes to test filtering)
    n_map_genes <- 40  # Less than total genes to test filtering
    phylomap <- data.frame(
        Stratum = sample(1:3, n_map_genes, replace = TRUE),
        GeneID = paste0("TestGene-", 1:n_map_genes)
    )
    
    # Test the function
    sc_set <- match_map_sc_seurat(
        seurat = seurat_obj,
        phylomap = phylomap,
        layer = "counts"
    )
    
    # Test that object was created correctly
    expect_s7_class(sc_set, myTAI::ScPhyloExpressionSet)
    expect_equal(sc_set@num_genes, n_map_genes)  # Should only include mapped genes
    expect_equal(sc_set@num_samples, n_cells)
    
    # Test that only genes in phylomap are retained
    expect_true(all(sc_set@gene_ids %in% phylomap$GeneID))
    
    # Test that strata were assigned correctly
    expect_equal(length(sc_set@strata), n_map_genes)
    expect_true(all(as.numeric(sc_set@strata) %in% phylomap$Stratum))
})
