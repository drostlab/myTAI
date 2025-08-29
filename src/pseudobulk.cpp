#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppThread.h>
// [[Rcpp::depends(RcppThread)]]
#include <thread>

//' @title Create Pseudobulk Expression Matrix (C++ Implementation)
//' @description Efficiently aggregate single-cell expression data by cell groups
//' to create pseudobulk expression profiles using optimized sparse matrix operations.
//' 
//' @param expression Sparse expression matrix (genes x cells) - dgCMatrix format
//' @param groups Integer vector indicating group membership for each cell (0-based indexing)
//' @param n_groups Integer, total number of unique groups
//' @param ncores Integer, number of cores to use for parallel processing (default: 10)
//' @return Dense matrix of pseudobulked expression (genes x groups)
//' 
//' @details
//' This function efficiently aggregates single-cell data by summing expression values
//' within each group. It processes genes in parallel and uses sparse matrix optimizations
//' to skip zero entries. The output is a dense matrix for consistency with downstream
//' analysis functions.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_pseudobulk(const arma::sp_mat& expression,
                         const arma::uvec& groups,
                         int n_groups,
                         int ncores = 10) {
    
    int n_genes = expression.n_rows;
    int n_cells = expression.n_cols;
    
    // Input validation
    if (groups.n_elem != n_cells) {
        Rcpp::stop("Length of groups must match number of cells in expression matrix");
    }
    
    if (arma::max(groups) >= n_groups) {
        Rcpp::stop("groups contains indices >= n_groups (remember: 0-based indexing)");
    }
    
    // Cap ncores at available cores
    int max_cores = std::thread::hardware_concurrency();
    if (max_cores == 0) max_cores = 1; // Fallback if detection fails
    
    if (ncores > max_cores) {
        Rcpp::warning("Requested %d cores but only %d cores detected. Using %d cores.", 
                      ncores, max_cores, max_cores);
        ncores = max_cores;
    }
    
    // Pre-allocate result matrix (genes x groups)
    arma::mat result(n_genes, n_groups, arma::fill::zeros);
    
    // Progress bar for genes
    RcppThread::ProgressBar bar(n_genes, 1);
    
    // Process genes in parallel
    RcppThread::parallelFor(0, n_genes, [&](int i) {
        // Get sparse row for this gene
        arma::sp_mat::const_row_iterator row_it = expression.begin_row(i);
        arma::sp_mat::const_row_iterator row_end = expression.end_row(i);
        
        // Aggregate expression values by group
        for (; row_it != row_end; ++row_it) {
            int cell_idx = row_it.col();
            double expr_val = *row_it;
            int group_idx = groups(cell_idx);
            result(i, group_idx) += expr_val;
        }
        
        // Update progress bar
        bar++;
    }, ncores);
    
    return result;
}

//' @title Create Pseudobulk Expression Matrix (Memory-Optimized Implementation)
//' @description Memory-optimized version that processes the pseudobulking in batches
//' for very large datasets to avoid memory issues.
//' 
//' @param expression Sparse expression matrix (genes x cells) - dgCMatrix format
//' @param groups Integer vector indicating group membership for each cell (0-based indexing)
//' @param n_groups Integer, total number of unique groups
//' @param batch_size Integer, number of genes to process per batch (default: 1000)
//' @param ncores Integer, number of cores to use for parallel processing (default: 10)
//' @return Dense matrix of pseudobulked expression (genes x groups)
//' 
//' @details
//' This function processes genes in batches to manage memory usage for very large
//' single-cell datasets. It's particularly useful when dealing with hundreds of
//' thousands of cells and tens of thousands of genes.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat cpp_pseudobulk_batched(const arma::sp_mat& expression,
                                 const arma::uvec& groups,
                                 int n_groups,
                                 int batch_size = 1000,
                                 int ncores = 10) {
    
    int n_genes = expression.n_rows;
    int n_cells = expression.n_cols;
    
    // Input validation
    if (groups.n_elem != n_cells) {
        Rcpp::stop("Length of groups must match number of cells in expression matrix");
    }
    
    if (arma::max(groups) >= n_groups) {
        Rcpp::stop("groups contains indices >= n_groups (remember: 0-based indexing)");
    }
    
    // Cap ncores at available cores
    int max_cores = std::thread::hardware_concurrency();
    if (max_cores == 0) max_cores = 1; // Fallback if detection fails
    
    if (ncores > max_cores) {
        Rcpp::warning("Requested %d cores but only %d cores detected. Using %d cores.", 
                      ncores, max_cores, max_cores);
        ncores = max_cores;
    }
    
    // Pre-allocate result matrix (genes x groups)
    arma::mat result(n_genes, n_groups, arma::fill::zeros);
    
    // Calculate number of batches
    int n_batches = (n_genes + batch_size - 1) / batch_size;
    
    // Progress bar for batches
    RcppThread::ProgressBar bar(n_batches, 1);
    
    // Process gene batches
    for (int batch_idx = 0; batch_idx < n_batches; batch_idx++) {
        // Calculate batch boundaries
        int start_gene = batch_idx * batch_size;
        int end_gene = std::min(start_gene + batch_size, n_genes);
        int batch_n_genes = end_gene - start_gene;
        
        // Extract batch submatrix
        arma::sp_mat batch_expr = expression.rows(start_gene, end_gene - 1);
        
        // Process genes in this batch in parallel
        RcppThread::parallelFor(0, batch_n_genes, [&](int local_i) {
            int global_i = start_gene + local_i;
            
            // Get sparse row for this gene
            arma::sp_mat::const_row_iterator row_it = batch_expr.begin_row(local_i);
            arma::sp_mat::const_row_iterator row_end = batch_expr.end_row(local_i);
            
            // Aggregate expression values by group
            for (; row_it != row_end; ++row_it) {
                int cell_idx = row_it.col();
                double expr_val = *row_it;
                int group_idx = groups(cell_idx);
                result(global_i, group_idx) += expr_val;
            }
        }, ncores);
        
        // Update progress bar
        bar++;
    }
    
    return result;
}
