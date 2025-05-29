
#include <random>
#include <iostream>
#include <algorithm>

#include <RcppParallel.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Core>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace RcppParallel;

struct PermuteWorker : public Worker {
    std::vector<double>& vec;
    Eigen::MatrixXd& mat;
    std::mt19937 g;
    
    PermuteWorker(std::vector<double>& vec_, Eigen::MatrixXd& mat_)
        : vec(vec_), mat(mat_), g(std::random_device{}()) {}
    
    void operator()(std::size_t begin, std::size_t end) {
        std::vector<double> local_vec = vec;
        
        for (std::size_t i = begin; i < end; i++) {
            std::shuffle(local_vec.begin(), local_vec.end(), g);
            Eigen::Map<Eigen::VectorXd> tmp_vec(local_vec.data(), local_vec.size());
            mat.col(i) = tmp_vec;
        }
    }
};

Eigen::MatrixXd permute_mat_parallel(Eigen::VectorXd vector, int num_permutations) {
    std::vector<double> vec(vector.data(), vector.data() + vector.size());
    Eigen::MatrixXd mat(vec.size(), num_permutations);
    PermuteWorker worker(vec, mat);
    parallelFor(0, mat.cols(), worker);
    return mat;
}


// [[Rcpp::export]]
Eigen::MatrixXd cpp_nullTXIs(Eigen::MatrixXd count_matrix, Eigen::VectorXd strata_vector, int num_permutations){
    // Normalise expression sample-wise
    Eigen::MatrixXd matrix_norm = count_matrix.array().rowwise() / count_matrix.colwise().sum().array();
    
    
    Eigen::MatrixXd perm_mat = permute_mat_parallel(strata_vector, num_permutations);
    Eigen::MatrixXd res = matrix_norm.transpose() * perm_mat;
    return res.transpose();
}







