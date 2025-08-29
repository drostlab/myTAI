# Test single-cell TXI calculation performance

test_that("C++ TXI implementation produces correct results", {
    skip_if_not_installed("Matrix")
    skip_on_cran()
    
    # Create a small test matrix for correctness testing
    set.seed(123)
    n_genes <- 100
    n_cells <- 50
    
    # Create sparse matrix with some zeros
    expr_data <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
    expr_data <- abs(expr_data) * 100  # Make values positive
    colnames(expr_data) <- paste0("cell_", 1:n_cells)
    rownames(expr_data) <- paste0("gene_", 1:n_genes)
    
    # Create strata values
    strata_values <- runif(n_genes, min = 1, max = 10)
    
    # R implementation (simplified version of original)
    txi_r <- function(expression, strata_values) {
        col_sums <- Matrix::colSums(expression)
        zero_cols <- col_sums == 0
        
        txi <- rep(NA_real_, ncol(expression))
        
        if (any(!zero_cols)) {
            non_zero_expr <- expression[, !zero_cols, drop = FALSE]
            non_zero_sums <- col_sums[!zero_cols]
            
            txi_non_zero <- as.numeric((Matrix::t(non_zero_expr) %*% strata_values) / non_zero_sums)
            txi[!zero_cols] <- txi_non_zero
        }
        
        names(txi) <- colnames(expression)
        return(txi)
    }
    
    # Calculate TXI with both methods
    txi_r_result <- txi_r(expr_data, strata_values)
    txi_cpp_result <- cpp_txi_sc(expr_data, strata_values, batch_size = 25, ncores = 1)
    
    # Also test the original .TXI_sc function
    txi_original_result <- .TXI_sc(expr_data, strata_values)
    
    # Test that results are nearly identical (allowing for floating point differences)
    expect_equal(length(txi_r_result), length(txi_cpp_result))
    expect_equal(length(txi_r_result), length(txi_original_result))
    
    # Check non-NA values are close
    non_na_r <- !is.na(unname(txi_r_result))
    non_na_cpp <- !is.na(unname(txi_cpp_result))
    non_na_original <- !is.na(unname(txi_original_result))
    expect_equal(non_na_r, non_na_cpp)
    expect_equal(non_na_r, non_na_original)
    
    if (any(non_na_r)) {
        expect_equal(unname(txi_r_result)[non_na_r], unname(txi_cpp_result)[non_na_cpp], tolerance = 1e-10)
        expect_equal(unname(txi_r_result)[non_na_r], unname(txi_original_result)[non_na_original], tolerance = 1e-10)
    }
    
})

test_that("TXI performance comparison (skipped by default)", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("microbenchmark")
    # skip("Performance test - run manually with testthat::test_file() and profile=TRUE")  # COMMENTED OUT TO RUN
    
    # Create larger test matrices for performance testing
    set.seed(42)
    
    test_sizes <- list(
        small = list(genes = 1000, cells = 500),
        medium = list(genes = 5000, cells = 2000),
        large = list(genes = 10000, cells = 5000),
        very_large = list(genes = 20000, cells = 100000),
        extreme = list(genes = 20000, cells = 200000)
    )
    
    results <- list()
    
    for (size_name in names(test_sizes)) {
        cat("\n=== Testing", size_name, "dataset ===\n")
        
        n_genes <- test_sizes[[size_name]]$genes
        n_cells <- test_sizes[[size_name]]$cells
        
        cat("Creating", n_genes, "genes x", n_cells, "cells sparse matrix...\n")
        
        # Create realistic sparse single-cell data
        expr_data <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.05)
        expr_data <- abs(expr_data) * rpois(length(expr_data@x), lambda = 5)
        colnames(expr_data) <- paste0("cell_", 1:n_cells)
        rownames(expr_data) <- paste0("gene_", 1:n_genes)
        
        # Create strata values
        strata_values <- sample(1:12, n_genes, replace = TRUE)
        
        cat("Matrix density:", round(length(expr_data@x) / (n_genes * n_cells) * 100, 2), "%\n")
        cat("Matrix size:", round(object.size(expr_data) / 1024^2, 1), "MB\n")
        
        # R implementation
        txi_r <- function(expression, strata_values) {
            col_sums <- Matrix::colSums(expression)
            zero_cols <- col_sums == 0
            
            txi <- rep(NA_real_, ncol(expression))
            
            if (any(!zero_cols)) {
                non_zero_expr <- expression[, !zero_cols, drop = FALSE]
                non_zero_sums <- col_sums[!zero_cols]
                
                txi_non_zero <- as.numeric((Matrix::t(non_zero_expr) %*% strata_values) / non_zero_sums)
                txi[!zero_cols] <- txi_non_zero
            }
            
            names(txi) <- colnames(expression)
            return(txi)
        }
        
        # For very large datasets, reduce the number of benchmark runs and focus on key implementations
        times <- if (n_cells > 150000) 2 else if (n_cells > 50000) 3 else 5
        
        # Skip some implementations for very large datasets to save time
        if (n_cells > 150000) {
            cat("Running minimal benchmarks for extreme dataset (200k+ cells)...\n")
            
            # Only test most promising configurations for extreme datasets
            bench_results <- microbenchmark::microbenchmark(
                original_TXI_sc = .TXI_sc(expr_data, strata_values),
                cpp_batched_8cores = cpp_txi_sc(expr_data, strata_values, batch_size = 2000, ncores = 8),
                cpp_batched_12cores = cpp_txi_sc(expr_data, strata_values, batch_size = 2000, ncores = 12),
                times = times,
                unit = "s"  # Use seconds for extreme datasets
            )
        } else if (n_cells > 50000) {
            cat("Running reduced benchmarks for very large dataset (100k cells) with high parallelism...\n")
            
            # Test with various core counts to find optimal parallelization
            bench_results <- microbenchmark::microbenchmark(
                original_TXI_sc = .TXI_sc(expr_data, strata_values),
                cpp_batched_4cores = cpp_txi_sc(expr_data, strata_values, batch_size = 2000, ncores = 4),
                cpp_batched_8cores = cpp_txi_sc(expr_data, strata_values, batch_size = 2000, ncores = 8),
                cpp_batched_12cores = cpp_txi_sc(expr_data, strata_values, batch_size = 2000, ncores = 12),
                times = times,
                unit = "s"  # Use seconds for large datasets
            )
        } else {
            # Also test higher core counts for smaller datasets to see scaling
            bench_results <- microbenchmark::microbenchmark(
                R_implementation = txi_r(expr_data, strata_values),
                original_TXI_sc = .TXI_sc(expr_data, strata_values),
                cpp_batched_1core = cpp_txi_sc(expr_data, strata_values, batch_size = 1000, ncores = 1),
                cpp_batched_4cores = cpp_txi_sc(expr_data, strata_values, batch_size = 1000, ncores = 4),
                cpp_batched_8cores = cpp_txi_sc(expr_data, strata_values, batch_size = 1000, ncores = 8),
                times = times,
                unit = "ms"
            )
        }
        
        print(bench_results)
        
        # Store results
        results[[size_name]] <- list(
            dimensions = c(genes = n_genes, cells = n_cells),
            density = length(expr_data@x) / (n_genes * n_cells),
            size_mb = as.numeric(object.size(expr_data) / 1024^2),
            benchmark = bench_results
        )
        
        # Test correctness
        cat("Verifying correctness...\n")
        txi_r_result <- txi_r(expr_data, strata_values)
        txi_original_result <- .TXI_sc(expr_data, strata_values)
        txi_cpp_result <- cpp_txi_sc(expr_data, strata_values, batch_size = 1000, ncores = 1)
        
        expect_equal(length(txi_r_result), length(txi_original_result))
        expect_equal(length(txi_r_result), length(txi_cpp_result))
        
        non_na_r <- !is.na(unname(txi_r_result))
        non_na_original <- !is.na(unname(txi_original_result))
        non_na_cpp <- !is.na(unname(txi_cpp_result))
        expect_equal(non_na_r, non_na_original)
        expect_equal(non_na_r, non_na_cpp)
        
        if (any(non_na_r)) {
            max_diff_original <- max(abs(unname(txi_r_result)[non_na_r] - unname(txi_original_result)[non_na_original]), na.rm = TRUE)
            max_diff_cpp <- max(abs(unname(txi_r_result)[non_na_r] - unname(txi_cpp_result)[non_na_cpp]), na.rm = TRUE)
            cat("Maximum difference between R and original .TXI_sc:", max_diff_original, "\n")
            cat("Maximum difference between R and C++ results:", max_diff_cpp, "\n")
            expect_true(max_diff_original < 1e-10, 
                       info = paste("Original results differ by more than 1e-10:", max_diff_original))
            expect_true(max_diff_cpp < 1e-10, 
                       info = paste("C++ results differ by more than 1e-10:", max_diff_cpp))
        }
        
        cat("Correctness verified!\n\n")
    }
    
    # Summary
    cat("\n=== PERFORMANCE SUMMARY ===\n")
    for (size_name in names(results)) {
        res <- results[[size_name]]
        bench <- res$benchmark
        
        cat("\n", toupper(size_name), "dataset (", 
            res$dimensions["genes"], " genes x ", 
            res$dimensions["cells"], " cells, ",
            round(res$density * 100, 2), "% density, ",
            round(res$size_mb, 1), " MB):\n", sep = "")
        
        medians <- aggregate(bench$time, by = list(bench$expr), median)
        names(medians) <- c("method", "median_time_ns")
        medians$median_time_ms <- medians$median_time_ns / 1e6
        
        # Use original .TXI_sc as baseline since it's the actual implementation
        baseline_time <- medians$median_time_ms[medians$method == "original_TXI_sc"]
        if (length(baseline_time) == 0) {
            # Fallback to R implementation if original not found
            baseline_time <- medians$median_time_ms[medians$method == "R_implementation"]
        }
        
        for (i in 1:nrow(medians)) {
            method <- medians$method[i]
            time_ms <- medians$median_time_ms[i]
            speedup <- baseline_time / time_ms
            
            baseline_indicator <- if (method == "original_TXI_sc") " (baseline)" else ""
            
            cat("  ", method, baseline_indicator, ": ", 
                round(time_ms, 1), " ms (", 
                round(speedup, 1), "x speedup)\n", sep = "")
        }
    }
    
    # Save results for further analysis
    saveRDS(results, file.path(tempdir(), "txi_benchmark_results.rds"))
    cat("\nBenchmark results saved to:", file.path(tempdir(), "txi_benchmark_results.rds"), "\n")
})

test_that("TXI handles edge cases correctly", {
    skip_if_not_installed("Matrix")
    skip_on_cran()
    
    # Test with all-zero matrix
    expr_zeros <- Matrix::Matrix(0, nrow = 10, ncol = 5, sparse = TRUE)
    strata_values <- 1:10
    
    txi_zeros <- cpp_txi_sc(expr_zeros, strata_values)
    expect_true(all(is.na(txi_zeros)))
    expect_equal(length(txi_zeros), 5)
    
    # Test with single non-zero cell
    expr_single <- Matrix::Matrix(0, nrow = 10, ncol = 5, sparse = TRUE)
    expr_single[1:5, 1] <- c(10, 20, 30, 40, 50)
    
    txi_single <- cpp_txi_sc(expr_single, strata_values)
    expect_equal(sum(!is.na(txi_single)), 1)  # Only one non-NA value
    
    expected_txi <- sum(c(10, 20, 30, 40, 50) * strata_values[1:5]) / sum(c(10, 20, 30, 40, 50))
    expect_equal(txi_single[1], expected_txi, tolerance = 1e-10)
    
    # Test with very large batch size
    expr_small <- Matrix::rsparsematrix(50, 20, density = 0.1)
    expr_small <- abs(expr_small) * 10
    strata_values_small <- runif(50, 1, 5)
    
    txi_large_batch <- cpp_txi_sc(expr_small, strata_values_small, batch_size = 1000)
    txi_small_batch <- cpp_txi_sc(expr_small, strata_values_small, batch_size = 5)
    
    expect_equal(txi_large_batch, txi_small_batch, tolerance = 1e-12)
})

test_that("Memory usage profiling (manual test)", {
    skip("Memory profiling test - run manually")
    skip_if_not_installed("profmem")
    skip_if_not_installed("Matrix")
    
    # Create medium-sized test data
    set.seed(123)
    n_genes <- 5000
    n_cells <- 2000
    
    expr_data <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.05)
    expr_data <- abs(expr_data) * rpois(length(expr_data@x), lambda = 3)
    strata_values <- sample(1:12, n_genes, replace = TRUE)
    
    cat("Testing memory usage for", n_genes, "genes x", n_cells, "cells\n")
    
    # Profile memory usage of C++ implementation
    mem_profile <- profmem::profmem({
        result <- cpp_txi_sc(expr_data, strata_values, batch_size = 500, ncores = 1)
    })
    
    cat("Memory profile:\n")
    print(mem_profile)
    
    total_mem <- sum(mem_profile$bytes, na.rm = TRUE)
    cat("Total memory allocated:", round(total_mem / 1024^2, 2), "MB\n")
    
    # Compare with different batch sizes
    cat("\nTesting different batch sizes:\n")
    batch_sizes <- c(100, 500, 1000, 2000, 5000)
    
    for (bs in batch_sizes) {
        if (bs <= n_cells) {
            mem_prof <- profmem::profmem({
                cpp_txi_sc(expr_data, strata_values, batch_size = bs, ncores = 1)
            })
            total <- sum(mem_prof$bytes, na.rm = TRUE)
            cat("Batch size", bs, ":", round(total / 1024^2, 2), "MB\n")
        }
    }
})
