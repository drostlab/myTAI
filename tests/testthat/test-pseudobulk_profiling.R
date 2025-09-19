# # Test pseudobulking performance with different implementations and dataset sizes

# library(testthat)
# library(Matrix)

# context("Pseudobulking Performance Profiling")

# # Test that C++ functions are available
# test_that("C++ pseudobulk functions are available", {
#     # Check if the functions exist
#     expect_true(exists("cpp_pseudobulk"))
#     expect_true(exists("cpp_pseudobulk_batched"))
    
#     # Test with minimal data
#     small_expr <- Matrix::sparseMatrix(i = c(1, 2), j = c(1, 1), x = c(5, 3), dims = c(2, 1))
#     small_groups <- c(0L)  # 0-based indexing
    
#     result1 <- cpp_pseudobulk(small_expr, small_groups, 1)
#     result2 <- cpp_pseudobulk_batched(small_expr, small_groups, 1, 1000)
    
#     expect_equal(dim(result1), c(2, 1))
#     expect_equal(dim(result2), c(2, 1))
# })

# test_that("Pseudobulk implementations produce correct results", {
#     skip_if_not_installed("Matrix")
#     skip_on_cran()
    
#     set.seed(123)
#     n_genes <- 100
#     n_cells <- 200
#     n_groups <- 5
    
#     # Create test data
#     expr_data <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.1)
#     expr_data <- abs(expr_data) * 10
#     groups <- factor(sample(1:n_groups, n_cells, replace = TRUE))
    
#     # R implementation
#     pseudobulk_r <- function(expression, groups) {
#         unique_groups <- levels(groups)
#         result_list <- lapply(unique_groups, function(group) {
#             group_cells <- which(groups == group)
#             if (length(group_cells) == 1) {
#                 expression[, group_cells, drop = FALSE]
#             } else {
#                 Matrix::rowSums(expression[, group_cells, drop = FALSE])
#             }
#         })
#         result <- do.call(cbind, result_list)
#         colnames(result) <- unique_groups
#         return(as.matrix(result))
#     }
    
#     # Test implementations
#     result_r <- pseudobulk_r(expr_data, groups)
    
#     # Convert groups to 0-based for C++
#     group_indices <- as.numeric(groups) - 1
#     result_cpp <- cpp_pseudobulk(expr_data, group_indices, n_groups)
#     result_cpp_batched <- cpp_pseudobulk_batched(expr_data, group_indices, n_groups, 50)
    
#     # Add column names to C++ results
#     colnames(result_cpp) <- levels(groups)
#     colnames(result_cpp_batched) <- levels(groups)
    
#     # Check correctness
#     expect_equal(dim(result_r), dim(result_cpp))
#     expect_equal(dim(result_r), dim(result_cpp_batched))
#     expect_equal(result_r, result_cpp, tolerance = 1e-10)
#     expect_equal(result_r, result_cpp_batched, tolerance = 1e-10)
# })

# test_that("Pseudobulk performance comparison", {
#     skip_if_not_installed("Matrix")
#     skip_if_not_installed("microbenchmark")
    
#     # Test scenarios similar to TXI profiling
#     test_sizes <- list(
#         small = list(genes = 1000, cells = 1000, groups = 5),
#         medium = list(genes = 5000, cells = 5000, groups = 8),
#         large = list(genes = 10000, cells = 20000, groups = 10),
#         very_large = list(genes = 20000, cells = 50000, groups = 15)
#     )
    
#     results <- list()
    
#     for (size_name in names(test_sizes)) {
#         cat("\n=== Testing", size_name, "dataset ===\n")
        
#         n_genes <- test_sizes[[size_name]]$genes
#         n_cells <- test_sizes[[size_name]]$cells
#         n_groups <- test_sizes[[size_name]]$groups
        
#         cat("Creating", n_genes, "genes x", n_cells, "cells sparse matrix...\n")
        
#         # Create test data
#         set.seed(42)
#         expr_data <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.05)
#         expr_data <- abs(expr_data) * rpois(length(expr_data@x), lambda = 5)
#         groups <- factor(sample(1:n_groups, n_cells, replace = TRUE))
        
#         cat("Matrix density:", round(length(expr_data@x) / (n_genes * n_cells) * 100, 2), "%\n")
#         cat("Matrix size:", round(object.size(expr_data) / 1024^2, 1), "MB\n")
        
#         # R implementation
#         pseudobulk_r <- function(expression, groups) {
#             unique_groups <- levels(groups)
#             result_list <- lapply(unique_groups, function(group) {
#                 group_cells <- which(groups == group)
#                 if (length(group_cells) == 1) {
#                     expression[, group_cells, drop = FALSE]
#                 } else {
#                     Matrix::rowSums(expression[, group_cells, drop = FALSE])
#                 }
#             })
#             result <- do.call(cbind, result_list)
#             colnames(result) <- unique_groups
#             return(as.matrix(result))
#         }
        
#         # Convert groups for C++
#         group_indices <- as.numeric(groups) - 1
        
#         # Benchmark runs
#         times <- if (n_cells > 30000) 3 else 5
        
#         bench_results <- microbenchmark::microbenchmark(
#             r_implementation = pseudobulk_r(expr_data, groups),
#             cpp_standard = cpp_pseudobulk(expr_data, group_indices, n_groups),
#             cpp_batched = cpp_pseudobulk_batched(expr_data, group_indices, n_groups, 1000),
#             times = times,
#             unit = if (n_cells > 10000) "s" else "ms"
#         )
        
#         print(bench_results)
        
#         # Store results
#         results[[size_name]] <- list(
#             dimensions = c(genes = n_genes, cells = n_cells, groups = n_groups),
#             density = length(expr_data@x) / (n_genes * n_cells),
#             size_mb = as.numeric(object.size(expr_data) / 1024^2),
#             benchmark = bench_results
#         )
        
#         # Verify correctness
#         cat("Verifying correctness...\n")
#         result_r <- pseudobulk_r(expr_data, groups)
#         result_cpp <- cpp_pseudobulk(expr_data, group_indices, n_groups)
#         result_cpp_batched <- cpp_pseudobulk_batched(expr_data, group_indices, n_groups, 1000)
        
#         expect_equal(dim(result_r), dim(result_cpp))
#         expect_equal(dim(result_r), dim(result_cpp_batched))
#         expect_equal(result_r, result_cpp, tolerance = 1e-10)
#         expect_equal(result_r, result_cpp_batched, tolerance = 1e-10)
        
#         cat("Correctness verified!\n\n")
#     }
    
#     # Summary
#     cat("\n=== PSEUDOBULKING PERFORMANCE SUMMARY ===\n\n")
#     for (size_name in names(results)) {
#         res <- results[[size_name]]
#         bench <- res$benchmark
        
#         cat(toupper(size_name), "dataset (", 
#             res$dimensions["genes"], " genes x ", 
#             res$dimensions["cells"], " cells x ",
#             res$dimensions["groups"], " groups, ",
#             round(res$density * 100, 2), "% density, ",
#             round(res$size_mb, 1), " MB):\n", sep = "")
        
#         medians <- aggregate(bench$time, by = list(bench$expr), median)
#         names(medians) <- c("method", "median_time_ns")
        
#         # Use R as baseline
#         baseline_time <- medians$median_time_ns[medians$method == "r_implementation"]
        
#         for (i in 1:nrow(medians)) {
#             method <- medians$method[i]
#             time_ns <- medians$median_time_ns[i]
#             speedup <- baseline_time / time_ns
            
#             # Convert to appropriate unit
#             if (max(medians$median_time_ns) > 1e9) {
#                 time_display <- paste(round(time_ns / 1e6, 1), "ms")
#             } else {
#                 time_display <- paste(round(time_ns / 1e6, 2), "ms")
#             }
            
#             baseline_indicator <- if (method == "r_implementation") " (baseline)" else ""
            
#             cat("  ", i, baseline_indicator, ": ", time_display, " (", 
#                 round(speedup, 1), "x speedup)\n", sep = "")
#         }
#         cat("\n")
#     }
    
#     # Save results
#     saveRDS(results, file.path(tempdir(), "pseudobulk_benchmark_results.rds"))
#     cat("Benchmark results saved to:", file.path(tempdir(), "pseudobulk_benchmark_results.rds"), "\n")
# })
