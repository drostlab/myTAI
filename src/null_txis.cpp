#include <RcppArmadillo/Lighter>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppThread.h>



//credit to Filipa Martins Costa

// [[Rcpp::export]]
arma::mat cpp_nullTXIs(const arma::mat& count_matrix, const arma::vec& strata_vector, int num_permutations) {
    int m = count_matrix.n_cols;

    arma::mat txis(num_permutations, m, arma::fill::zeros);
    
    RcppThread::ProgressBar bar(num_permutations, 1);
    
    arma::mat count_norm = count_matrix.each_row() / sum(count_matrix, 0);
    
    RcppThread::parallelFor(0, num_permutations, [&](int i) {
        arma::vec perm = shuffle(strata_vector);
        txis.row(i) = sum(count_norm.each_col() % perm, 0);
        bar++;
    });
    
    
    return txis;
}
