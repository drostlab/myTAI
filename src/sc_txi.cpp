#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppThread.h>
// [[Rcpp::depends(RcppThread)]]
#include <thread>

//' @title Calculate TXI for Single-Cell Expression Data (C++ Implementation)
//' @description Efficiently calculate TXI values for sparse single-cell expression matrices
//' using batch processing and parallel computation.
//' 
//' @param expression Sparse expression matrix (genes x cells) - dgCMatrix format
//' @param strata_values Numeric vector of phylostratum values for each gene
//' @param batch_size Integer, number of cells to process per batch (default: 2000)
//' @param ncores Integer, number of cores to use for parallel processing (default: 10, automatically capped at available cores)
//' @return Numeric vector of TXI values for each cell
//' 
//' @details
//' This function processes large sparse single-cell expression matrices efficiently by:
//' - Splitting cells into batches to manage memory usage
//' - Using parallel processing across batches
//' - Leveraging sparse matrix operations to skip zero entries
//' - Handling cells with zero expression by returning NA
//' 
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_txi_sc(const arma::sp_mat& expression, 
                               const arma::vec& strata_values,
                               int batch_size = 2000,
                               int ncores = 10) {
    
    int n_cells = expression.n_cols;
    int n_genes = expression.n_rows;
    
    // Input validation
    if (strata_values.n_elem != n_genes) {
        Rcpp::stop("Length of strata_values must match number of genes in expression matrix");
    }
    
    // Cap ncores at available cores
    int max_cores = std::thread::hardware_concurrency();
    if (max_cores == 0) max_cores = 1; // Fallback if detection fails
    
    if (ncores > max_cores) {
        Rcpp::warning("Requested %d cores but only %d cores detected. Using %d cores.", 
                      ncores, max_cores, max_cores);
        ncores = max_cores;
    }
    
    // Pre-allocate result vector
    Rcpp::NumericVector txi(n_cells);
    std::fill(txi.begin(), txi.end(), NA_REAL);
    
    // Calculate number of batches
    int n_batches = (n_cells + batch_size - 1) / batch_size;
    
    // Progress bar for batches
    RcppThread::ProgressBar bar(n_batches, 1);
    
    // Process batches in parallel
    RcppThread::parallelFor(0, n_batches, [&](int batch_idx) {
        // Calculate batch boundaries
        int start_col = batch_idx * batch_size;
        int end_col = std::min(start_col + batch_size, n_cells);
        int batch_n_cols = end_col - start_col;
        
        // Extract batch submatrix
        arma::sp_mat batch_expr = expression.cols(start_col, end_col - 1);
        
        // Calculate column sums for this batch manually
        arma::rowvec col_sums(batch_n_cols, arma::fill::zeros);
        for (int j = 0; j < batch_n_cols; j++) {
            arma::sp_mat::const_col_iterator col_it = batch_expr.begin_col(j);
            arma::sp_mat::const_col_iterator col_end = batch_expr.end_col(j);
            double sum = 0.0;
            for (; col_it != col_end; ++col_it) {
                sum += *col_it;
            }
            col_sums(j) = sum;
        }
        
        // Process each cell in the batch
        for (int j = 0; j < batch_n_cols; j++) {
            double col_sum = col_sums(j);
            
            // Skip cells with zero expression
            if (col_sum == 0.0) {
                continue; // txi[start_col + j] remains NA
            }
            
            // Calculate weighted sum: sum(expression * strata_values)
            double weighted_sum = 0.0;
            
            // Iterate over non-zero elements in this column
            arma::sp_mat::const_col_iterator col_it = batch_expr.begin_col(j);
            arma::sp_mat::const_col_iterator col_end = batch_expr.end_col(j);
            
            for (; col_it != col_end; ++col_it) {
                int gene_idx = col_it.row();
                double expr_val = *col_it;
                weighted_sum += expr_val * strata_values(gene_idx);
            }
            
            // Calculate TXI
            txi[start_col + j] = weighted_sum / col_sum;
        }
        
        // Update progress bar
        bar++;
    }, ncores);
    
    // Set names if the expression matrix has column names
    // Note: Armadillo sparse matrices don't store colnames, 
    // so this would need to be handled at the R level
    
    return txi;
}