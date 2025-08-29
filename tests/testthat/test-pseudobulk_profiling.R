# Test pseudobulking performance with different implementations and dataset sizes

# Required packages
if (!require("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!require("microbenchmark", quietly = TRUE)) install.packages("microbenchmark")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(testthat)
library(Matrix)
library(microbenchmark)
library(ggplot2)

context("Pseudobulking Performance Profiling")

# Helper function to create synthetic single-cell data
create_synthetic_sc_data <- function(n_genes, n_cells, n_groups = 5, sparsity = 0.95, seed = 123) {
  set.seed(seed)
  
  # Create sparse expression matrix
  n_nonzero <- as.integer((1 - sparsity) * n_genes * n_cells)
  
  # Random positions for non-zero entries
  gene_indices <- sample(1:n_genes, n_nonzero, replace = TRUE)
  cell_indices <- sample(1:n_cells, n_nonzero, replace = TRUE)
  
  # Random expression values (log-normal-ish distribution)
  values <- rpois(n_nonzero, lambda = 3) + 1
  
  # Create sparse matrix
  expression <- sparseMatrix(
    i = gene_indices,
    j = cell_indices, 
    x = values,
    dims = c(n_genes, n_cells)
  )
  
  # Add row and column names
  rownames(expression) <- paste0("Gene_", 1:n_genes)
  colnames(expression) <- paste0("Cell_", 1:n_cells)
  
  # Create groups (factor)
  groups <- factor(sample(paste0("Group_", 1:n_groups), n_cells, replace = TRUE))
  
  return(list(expression = expression, groups = groups))
}

# Rebuild C++ functions after adding new file
test_that("C++ functions compile successfully", {
  expect_silent({
    # This will trigger recompilation if needed
    devtools::load_all()
  })
})

# Original R implementation for comparison
pseudobulk_r_original <- function(expression, groups) {
  # Use all unique identities present (levels of the factor)
  unique_groups <- levels(groups)
  
  # Pseudobulk by summing within each group
  pseudobulk_list <- lapply(unique_groups, function(group) {
    group_cells <- which(groups == group)
    if (length(group_cells) == 1) {
      expression[, group_cells, drop = FALSE]
    } else {
      if (inherits(expression, "sparseMatrix")) {
        Matrix::rowSums(expression[, group_cells, drop = FALSE])
      } else {
        rowSums(expression[, group_cells, drop = FALSE])
      }
    }
  })
  
  # If only one cell per type, ensure output is a matrix
  result <- do.call(cbind, pseudobulk_list)
  colnames(result) <- unique_groups
  
  # Convert to regular matrix for consistency with base class
  return(as.matrix(result))
}

# Optimized R implementation using Matrix operations
pseudobulk_r_optimized <- function(expression, groups) {
  unique_groups <- levels(groups)
  n_groups <- length(unique_groups)
  
  # Create group indicator matrix (sparse)
  group_indicators <- Matrix::sparseMatrix(
    i = as.numeric(groups),
    j = seq_along(groups),
    x = 1,
    dims = c(n_groups, ncol(expression))
  )
  rownames(group_indicators) <- unique_groups
  
  # Matrix multiplication for pseudobulking
  result <- as.matrix(expression %*% Matrix::t(group_indicators))
  colnames(result) <- unique_groups
  
  return(result)
}

# Wrapper functions for C++ implementations
cpp_pseudobulk_wrapper <- function(expression, groups) {
  # Convert groups to 0-based integer indices
  group_levels <- levels(groups)
  group_indices <- as.numeric(groups) - 1  # Convert to 0-based
  n_groups <- length(group_levels)
  
  result <- cpp_pseudobulk(expression, group_indices, n_groups)
  colnames(result) <- group_levels
  
  return(result)
}

cpp_pseudobulk_batched_wrapper <- function(expression, groups, batch_size = 1000) {
  group_levels <- levels(groups)
  group_indices <- as.numeric(groups) - 1
  n_groups <- length(group_levels)
  
  result <- cpp_pseudobulk_batched(expression, group_indices, n_groups, batch_size)
  colnames(result) <- group_levels
  
  return(result)
}

cpp_pseudobulk_cell_parallel_wrapper <- function(expression, groups) {
  group_levels <- levels(groups)
  group_indices <- as.numeric(groups) - 1
  n_groups <- length(group_levels)
  
  result <- cpp_pseudobulk_cell_parallel(expression, group_indices, n_groups)
  colnames(result) <- group_levels
  
  return(result)
}

# Test correctness on small dataset
test_that("All implementations produce identical results", {
  data <- create_synthetic_sc_data(n_genes = 100, n_cells = 500, n_groups = 3, sparsity = 0.8)
  
  result_r_orig <- pseudobulk_r_original(data$expression, data$groups)
  result_r_opt <- pseudobulk_r_optimized(data$expression, data$groups)
  result_cpp <- cpp_pseudobulk_wrapper(data$expression, data$groups)
  result_cpp_batched <- cpp_pseudobulk_batched_wrapper(data$expression, data$groups, batch_size = 50)
  result_cpp_cell <- cpp_pseudobulk_cell_parallel_wrapper(data$expression, data$groups)
  
  # Check dimensions
  expect_equal(dim(result_r_orig), dim(result_r_opt))
  expect_equal(dim(result_r_orig), dim(result_cpp))
  expect_equal(dim(result_r_orig), dim(result_cpp_batched))
  expect_equal(dim(result_r_orig), dim(result_cpp_cell))
  
  # Check values (allowing for small numerical differences)
  expect_equal(result_r_orig, result_r_opt, tolerance = 1e-10)
  expect_equal(result_r_orig, result_cpp, tolerance = 1e-10)
  expect_equal(result_r_orig, result_cpp_batched, tolerance = 1e-10)
  expect_equal(result_r_orig, result_cpp_cell, tolerance = 1e-10)
  
  # Check column names
  expect_equal(colnames(result_r_orig), colnames(result_cpp))
  expect_equal(colnames(result_r_orig), colnames(result_cpp_batched))
  expect_equal(colnames(result_r_orig), colnames(result_cpp_cell))
})

# Performance benchmarking with different dataset sizes
test_that("Performance comparison across dataset sizes", {
  
  # Define test scenarios
  scenarios <- list(
    small = list(n_genes = 1000, n_cells = 5000, n_groups = 5),
    medium = list(n_genes = 5000, n_cells = 20000, n_groups = 8),
    large = list(n_genes = 10000, n_cells = 50000, n_groups = 10),
    very_large = list(n_genes = 20000, n_cells = 100000, n_groups = 15),
    extreme = list(n_genes = 30000, n_cells = 200000, n_groups = 20)
  )
  
  results_list <- list()
  
  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    cat("\nTesting scenario:", scenario_name, "\n")
    cat("Genes:", scenario$n_genes, "Cells:", scenario$n_cells, "Groups:", scenario$n_groups, "\n")
    
    # Create test data
    data <- create_synthetic_sc_data(
      n_genes = scenario$n_genes, 
      n_cells = scenario$n_cells, 
      n_groups = scenario$n_groups,
      sparsity = 0.93  # High sparsity like real sc-RNA-seq
    )
    
    cat("Matrix sparsity:", (1 - length(data$expression@x) / (nrow(data$expression) * ncol(data$expression))), "\n")
    
    # For very large datasets, reduce the number of benchmark runs
    times <- if (scenario$n_cells > 50000) 3 else 5
    
    # Run benchmarks
    benchmark_result <- microbenchmark(
      r_original = pseudobulk_r_original(data$expression, data$groups),
      r_optimized = pseudobulk_r_optimized(data$expression, data$groups),
      cpp_standard = cpp_pseudobulk_wrapper(data$expression, data$groups),
      cpp_batched = cpp_pseudobulk_batched_wrapper(data$expression, data$groups, batch_size = 1000),
      cpp_cell_parallel = cpp_pseudobulk_cell_parallel_wrapper(data$expression, data$groups),
      times = times,
      unit = "s"
    )
    
    # Add scenario information
    benchmark_result$scenario <- scenario_name
    benchmark_result$n_genes <- scenario$n_genes
    benchmark_result$n_cells <- scenario$n_cells
    benchmark_result$n_groups <- scenario$n_groups
    
    results_list[[scenario_name]] <- benchmark_result
    
    # Print summary for this scenario
    print(summary(benchmark_result))
    cat("\n")
  }
  
  # Combine all results
  all_results <- do.call(rbind, results_list)
  
  # Create comprehensive visualization
  p1 <- ggplot(all_results, aes(x = scenario, y = time/1e9, fill = expr)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(
      title = "Pseudobulking Performance Comparison",
      subtitle = "Across Different Dataset Sizes",
      x = "Dataset Size", 
      y = "Time (seconds, log scale)",
      fill = "Implementation"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p1)
  
  # Create speedup comparison relative to original R implementation
  speedup_data <- all_results %>%
    group_by(scenario, n_genes, n_cells) %>%
    summarise(
      r_original_median = median(time[expr == "r_original"]),
      r_optimized_speedup = r_original_median / median(time[expr == "r_optimized"]),
      cpp_standard_speedup = r_original_median / median(time[expr == "cpp_standard"]),
      cpp_batched_speedup = r_original_median / median(time[expr == "cpp_batched"]),
      cpp_cell_parallel_speedup = r_original_median / median(time[expr == "cpp_cell_parallel"]),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = ends_with("_speedup"),
      names_to = "implementation",
      values_to = "speedup"
    ) %>%
    mutate(
      implementation = str_remove(implementation, "_speedup"),
      dataset_label = paste0(scenario, "\n(", scales::comma(n_genes), " genes, ", scales::comma(n_cells), " cells)")
    )
  
  p2 <- ggplot(speedup_data, aes(x = dataset_label, y = speedup, fill = implementation)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(
      title = "Speedup Relative to Original R Implementation",
      subtitle = "Values > 1 indicate faster performance",
      x = "Dataset Size", 
      y = "Speedup (x times faster)",
      fill = "Implementation"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(trans = "log10", labels = function(x) paste0(x, "x"))
  
  print(p2)
  
  # Memory usage comparison (approximate)
  memory_data <- all_results %>%
    group_by(scenario, expr, n_genes, n_cells) %>%
    summarise(median_time = median(time), .groups = "drop") %>%
    mutate(
      approx_memory_gb = (n_genes * n_cells * 8) / (1024^3),  # Approximate for dense matrix
      dataset_label = paste0(scenario, "\n(", scales::comma(n_genes), " genes, ", scales::comma(n_cells), " cells)")
    )
  
  p3 <- ggplot(memory_data, aes(x = approx_memory_gb, y = median_time/1e9, color = expr)) +
    geom_point(size = 3) +
    geom_line(aes(group = expr), alpha = 0.7) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Performance vs Dataset Size",
      subtitle = "Time complexity analysis",
      x = "Approximate Dense Matrix Size (GB, log scale)", 
      y = "Median Time (seconds, log scale)",
      color = "Implementation"
    ) +
    theme_minimal()
  
  print(p3)
  
  # Store results for later analysis
  assign("pseudobulk_benchmark_results", all_results, envir = .GlobalEnv)
  
  # Test that C++ implementations are consistently faster than R for large datasets
  large_scenario_results <- all_results[all_results$scenario %in% c("large", "very_large", "extreme"), ]
  
  if (nrow(large_scenario_results) > 0) {
    median_times <- large_scenario_results %>%
      group_by(scenario, expr) %>%
      summarise(median_time = median(time), .groups = "drop")
    
    for (scenario in c("large", "very_large", "extreme")) {
      scenario_times <- median_times[median_times$scenario == scenario, ]
      if (nrow(scenario_times) > 0) {
        r_orig_time <- scenario_times$median_time[scenario_times$expr == "r_original"]
        cpp_times <- scenario_times$median_time[grepl("cpp_", scenario_times$expr)]
        
        if (length(r_orig_time) > 0 && length(cpp_times) > 0) {
          expect_true(all(cpp_times < r_orig_time), 
                     info = paste("C++ implementations should be faster than R for scenario:", scenario))
        }
      }
    }
  }
})

# Edge case testing
test_that("Edge cases handled correctly", {
  
  # Single cell per group
  data_single <- create_synthetic_sc_data(n_genes = 100, n_cells = 5, n_groups = 5, sparsity = 0.7)
  result_single <- cpp_pseudobulk_wrapper(data_single$expression, data_single$groups)
  expect_equal(ncol(result_single), 5)
  expect_equal(nrow(result_single), 100)
  
  # Single group (all cells in one group)
  single_group <- factor(rep("Group_1", 1000))
  data_one_group <- create_synthetic_sc_data(n_genes = 500, n_cells = 1000, n_groups = 1, sparsity = 0.8)
  data_one_group$groups <- single_group
  result_one_group <- cpp_pseudobulk_wrapper(data_one_group$expression, data_one_group$groups)
  expect_equal(ncol(result_one_group), 1)
  expect_equal(colnames(result_one_group), "Group_1")
  
  # Very sparse matrix (99% zeros)
  data_sparse <- create_synthetic_sc_data(n_genes = 1000, n_cells = 5000, n_groups = 5, sparsity = 0.99)
  result_sparse <- cpp_pseudobulk_wrapper(data_sparse$expression, data_sparse$groups)
  expect_equal(dim(result_sparse), c(1000, 5))
})

# Final summary
test_that("Performance summary", {
  if (exists("pseudobulk_benchmark_results")) {
    cat("\n=== PSEUDOBULKING PERFORMANCE SUMMARY ===\n")
    
    summary_stats <- pseudobulk_benchmark_results %>%
      group_by(expr) %>%
      summarise(
        min_time_s = min(time) / 1e9,
        median_time_s = median(time) / 1e9,
        max_time_s = max(time) / 1e9,
        .groups = "drop"
      ) %>%
      arrange(median_time_s)
    
    print(summary_stats)
    
    # Calculate overall speedup
    r_orig_median <- summary_stats$median_time_s[summary_stats$expr == "r_original"]
    speedups <- r_orig_median / summary_stats$median_time_s
    names(speedups) <- summary_stats$expr
    
    cat("\nOverall speedup relative to original R implementation:\n")
    for (impl in names(speedups)) {
      cat(sprintf("%s: %.1fx faster\n", impl, speedups[[impl]]))
    }
    
    cat("\n=== RECOMMENDATIONS ===\n")
    cat("For datasets with:\n")
    cat("- < 10,000 cells: Use cpp_standard (best balance of speed and simplicity)\n")
    cat("- 10,000-50,000 cells: Use cpp_standard or cpp_cell_parallel\n") 
    cat("- > 50,000 cells: Use cpp_batched for memory efficiency\n")
    cat("- Very wide datasets (many genes, fewer cells): Use cpp_standard\n")
    cat("- Very tall datasets (fewer genes, many cells): Use cpp_cell_parallel\n")
  }
})
