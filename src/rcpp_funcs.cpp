//#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#include <math.h>
#include <map>
#include <random>
#include <iostream>
#include <RcppThread.h>
#include <string.h>

#ifdef _OPENMP
// OpenMP is available
// Include multiprocessing libraries and use parallelization
  #include <omp.h>
#endif
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Core>
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;
using namespace arma;



/* This whole 'permut' function has been adapted and taken from:
http://gallery.rcpp.org/articles/stl-random-shuffle/
- wrapper around R's RNG such that we get a uniform distribution over
- [0,n) as required by the STL algorithm
*/

int randWrapper(const int& n)
{ 
        return floor(unif_rand() * n); 
}

// Initializing the random number generator outside of the function
std::random_device rng;

std::default_random_engine gn(rng());
std::minstd_rand mrgn(42);
std::mt19937_64 urng(rng());


void updateProgressBar(int currentProgress, int totalProgress, int barWidth = 40) {
  float progressRatio = static_cast<float>(currentProgress) / totalProgress;
  int completedWidth = static_cast<int>(progressRatio * barWidth);
  
  std::string progressBar;
  progressBar.reserve(barWidth + 6);
  progressBar += '[';
  for (int i = 0; i < completedWidth; ++i) {
    progressBar += '=';
  }
  if (completedWidth < barWidth) {
    progressBar += '>';
    for (int i = completedWidth + 1; i < barWidth; ++i) {
      progressBar += ' ';
    }
  } else {
    progressBar += '=';
  }
  progressBar += ']';
  
  float percentage = round(progressRatio * 10000)/100;
  Rcout << '\r' << progressBar << " " << percentage << "%   ";
  // std::cout.flush();
}


Eigen::MatrixXd permut_mat(const Eigen::VectorXd& a,const int& permutations) {
  // already added by sourceCpp(), but needed standalone
  Eigen::MatrixXd permutedMat(permutations, a.size());
  Eigen::VectorXd shuffledVec = a;
  
  const int updateFrequency = 200;  // Update progress every x iterations
#ifdef _OPENMP  
  std::atomic<int> progress(0);
  
#pragma omp parallel
{
  int localProgress = 0;
  std::mt19937_64 urngp(rng());
  Eigen::VectorXd shuff = a;
  Eigen::VectorXd shuffl = shuff;

#pragma omp for
  for (int i = 0; i < permutations; i++) {
    shuffl = shuff;
    std::shuffle(shuffl.data(), shuffl.data() + shuffl.size(), urngp);
    permutedMat.row(i) = shuffl.transpose();
    localProgress++;
    
    
    if (localProgress % updateFrequency == 0) {
      progress.fetch_add(updateFrequency, std::memory_order_relaxed);
      
      // Display progress
      
      int currentProgress = progress.load(std::memory_order_relaxed);
      
      updateProgressBar(currentProgress,permutations);
    }
  }
}

#else
for (int i = 0; i < permutations; i++) {
  shuffledVec = a;
  std::shuffle(shuffledVec.data(), shuffledVec.data() + shuffledVec.size(), urng);
  permutedMat.row(i) = shuffledVec.transpose();
  if (i % updateFrequency == 0){
    updateProgressBar(i,permutations);
  }
  
}
#endif
updateProgressBar(permutations,permutations);
return permutedMat;
}

NumericVector permut(const NumericVector& a)
{
  
  // already added by sourceCpp(), but needed standalone
  // RNGScope scope;             
  
  // clone a into b to leave a alone
  NumericVector b = clone(a);
  
  std::shuffle(b.begin(), b.end(), urng);
  
  return b;
}

// @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_TAI(const Eigen::MatrixXd& ExpressionMatrix, const Eigen::VectorXd& Phylostratum) {
  
  Eigen::VectorXd Divisor = ExpressionMatrix.colwise().sum();
  Eigen::MatrixXd fMatrix = ExpressionMatrix.array().rowwise() / Divisor.transpose().array();
  Eigen::VectorXd total = Phylostratum.transpose() * fMatrix;
  
  return total;
}
/*
// @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_TAI_par(const Eigen::SparseMatrix<double> & ExpressionMatrix, const Eigen::VectorXd& Phylostratum) {
  
  Eigen::VectorXd Divisor = ExpressionMatrix.transpose() * Eigen::VectorXd::Ones(ExpressionMatrix.rows());
  Eigen::MatrixXd phylExp = Phylostratum.replicate(1, 7);
  phylExp = phylExp.array().rowwise() / Divisor.transpose().array();
  Eigen::VectorXd total = phylExp * ExpressionMatrix;
  
  return total;
}
*/

// @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_bootMatrix(const Eigen::MatrixXd& ExpressionMatrix, const Eigen::VectorXd& AgeVector, const int& permutations) 
{
        Rcout << std::endl;
        Rcout << "[ Number of Eigen threads that are employed on your machine: " << Eigen::nbThreads() << " ]" << std::endl;
        Rcout << std::endl;
        Eigen::VectorXd Divisor = ExpressionMatrix.colwise().sum();
        Eigen::MatrixXd fMatrix = ExpressionMatrix.array().rowwise() / Divisor.transpose().array();
        Rcout << "[ Computing age assignment permutations for test statistic ..." << " ]" << std::endl;
        Eigen::MatrixXd permMatrix = permut_mat(AgeVector,permutations);
        Rcout << std::endl;
        Rcout << "[ Computing variances of permuted transcriptome signatures ..." << " ]" << std::endl;
        Eigen::MatrixXd bootM = permMatrix * fMatrix;
        Rcout << std::endl;
        return bootM;
}

// @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_pMatrix(const Eigen::MatrixXd& ExpressionSet,const Eigen::VectorXd& AgeVector)
{
  Eigen::VectorXd Divisor = ExpressionSet.colwise().sum();
  Eigen::MatrixXd results = (ExpressionSet.array().rowwise() / Divisor.transpose().array()).colwise() * AgeVector.array();
  return results;
  
}


// @export
// [[Rcpp::export]]
double cpp_std_error(const NumericVector& x)
{
        
        return sd(x)/sqrt(x.size());
        
}

// @export 
// [[Rcpp::export]] 
double cpp_geom_mean(const NumericVector& x)
{
        
        return exp(mean(log(x)));
        
} 


// @export
// [[Rcpp::export]] 
double cpp_harmonic_mean(const NumericVector& x)
{
        double sum_val = 0.0;
        
        for (int i = 0; i < x.size(); i++){
                sum_val += 1.0/x[i];
        }
        return x.size() / sum_val;
} 



// @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_omitMatrix(const Eigen::MatrixXd& ExpressionSet, const Eigen::VectorXd& AgeVector){
  
  int nRows = ExpressionSet.rows();
  Eigen::VectorXd Divisor = ExpressionSet.colwise().sum();
  Eigen::VectorXd Numerator = AgeVector.transpose() * ExpressionSet;
  Eigen::MatrixXd AgeWeighted = ExpressionSet.array().colwise() * AgeVector.array();
  Eigen::MatrixXd Numer = Numerator.transpose().replicate(nRows, 1) - AgeWeighted;
  Eigen::MatrixXd ResultMatrix = Numer.cwiseQuotient(Divisor.transpose().replicate(nRows, 1) - ExpressionSet);
  return ResultMatrix;
}



//' @title rcpp_tei_parallel
//' @name rcpp_tei_parallel
//' @import Rcpp
//' @import Matrix
//' @description computes the phylogenetically based
//' transcriptome evolutionary index (TEI)
//' @return list
//' @param expression ExpressionSet as sparseMatrix
//' @param ps named Phylostratum
//' @param ncores number of cores
//' @examples
//' ## load example sequence data
//' data("PhyloExpressionSetExample", package="myTAI")
//' spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]), "sparseMatrix")
//' rownames(spmat) <- PhyloExpressionSetExample$GeneID
//' ps <- setNames(PhyloExpressionSetExample$Phylostratum, PhyloExpressionSetExample$GeneID)
//' rcpp_tei_parallel(spmat, ps)
//' @author Kristian K Ullrich
//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_tei_parallel(const arma::sp_mat& expression,
                             Rcpp::NumericVector ps,
                             int ncores = 1){
  std::vector< std::string > psnames =  ps.attr("names");
  int n_col = expression.n_cols;
  Rcpp::NumericVector sumx(n_col);
  Rcpp::NumericVector teisum(n_col);
  Rcpp::NumericVector tei(n_col);
  RcppThread::ProgressBar bar(n_col, 1);
  RcppThread::parallelFor(0, n_col, [&] (int j) {
    for (size_t i = 0; i < expression.n_rows; i++) {
      sumx[j] += expression(i, j);
      teisum[j] += expression(i, j) * ps[i];
      tei[j] = teisum[j]/sumx[j];
    }
  }, ncores);
  return Rcpp::List::create(Rcpp::Named("sumx") = sumx,
                            Rcpp::Named("teisum") = teisum, Rcpp::Named("tei") = tei);
}

//' @import Rcpp
//' @import Matrix
//' @title rcpp_boottei_parallel
//' @name rcpp_boottei_parallel
//' @description computes the phylogenetically based
//' transcriptome evolutionary index (TEI) shuffling the strata for permutation
//' statistic
//' @return sparseMatrix
//' @param expression ExpressionSet as sparseMatrix
//' @param ps named Phylostratum
//' @param permutations number of permutations
//' @param ncores number of cores
//' @examples
//' ## load example PhyloExpressionSetExample
//'
//' data("PhyloExpressionSetExample", package="myTAI")
//'
//' ## convert into sparseMatrix - rownames GeneID
//'
//' spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]),
//'     "sparseMatrix")
//' rownames(spmat) <- PhyloExpressionSetExample$GeneID
//'
//' ## create named Phylostratum vector
//'
//' ps <- setNames(PhyloExpressionSetExample$Phylostratum,
//'     PhyloExpressionSetExample$GeneID)
//'
//' ## get permutations
//' rcpp_boottei_parallel(spmat, ps, 100, 1)
//' @author Kristian K Ullrich
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_boottei_parallel(const arma::sp_mat& expression,
                                          Rcpp::NumericVector ps,
                                          const int& permutations,
                                          int ncores = 1){
  std::vector< std::string > psnames =  ps.attr("names");
  int n_col = expression.n_cols;
  int n_row = expression.n_rows;
  Rcpp::NumericMatrix fMatrix(n_row, n_col);
  Rcpp::NumericVector sampledVector(ps.length());
  Rcpp::NumericMatrix sampledMatrix(ps.length(), permutations);
  for (size_t l = 0; l < permutations; l++) {
    sampledVector = permut(ps);
    for(size_t m = 0; m < ps.length(); m++) {
      sampledMatrix(m, l) = sampledVector[m];
    }
  }
  Rcpp::NumericVector sumx(n_col);
  Rcpp::NumericMatrix bootM(permutations,n_col);
  RcppThread::ProgressBar bar(n_col, 1);
  RcppThread::parallelFor(0, n_col, [&] (int j) {
    for (size_t i = 0; i < n_row; i++) {
      sumx[j] += expression(i, j);
    }
    for (size_t k = 0; k < n_row; k++) {
      fMatrix(k, j) = expression(k, j) / sumx[j];
    }
    for (size_t n = 0; n < permutations; n++) {
      double teisum = 0;
      for (size_t o = 0; o < n_row; o++) {
        teisum += (sampledMatrix(o, n) * fMatrix(o, j));
      }
      bootM(n, j) = teisum;
    }
  }, ncores);
  return bootM;
}

//' @import Rcpp
//' @import Matrix
//' @title rcpp_pMatrix_parallel
//' @name rcpp_pMatrix_parallel
//' @description computes the partial
//' transcriptome evolutionary index (TEI) values for each single gene
//' @return sparseMatrix
//' @param expression ExpressionSet as sparseMatrix
//' @param ps named Phylostratum
//' @param ncores number of cores
//' @examples
//' ## load example PhyloExpressionSetExample
//'
//' data("PhyloExpressionSetExample", package="myTAI")
//'
//' ## convert into sparseMatrix - rownames GeneID
//'
//' spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]),
//'     "sparseMatrix")
//' rownames(spmat) <- PhyloExpressionSetExample$GeneID
//'
//' ## create named Phylostratum vector
//'
//' ps <- setNames(PhyloExpressionSetExample$Phylostratum,
//'     PhyloExpressionSetExample$GeneID)
//'
//' ## get pMatrix
//' rcpp_pMatrix_parallel(spmat, ps)
//' @author Kristian K Ullrich
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_pMatrix_parallel(const arma::sp_mat& expression,
                                          Rcpp::NumericVector ps,
                                          int ncores = 1){
  std::vector< std::string > psnames =  ps.attr("names");
  int n_col = expression.n_cols;
  int n_row = expression.n_rows;
  Rcpp::NumericMatrix pMatrix(n_row, n_col);
  Rcpp::NumericVector sumx(n_col);
  RcppThread::ProgressBar bar(n_col, 1);
  RcppThread::parallelFor(0, n_col, [&] (int j) {
    for (size_t i = 0; i < n_row; i++) {
      sumx[j] += expression(i, j);
    }
    for (size_t k = 0; k < n_row; k++) {
      pMatrix(k, j) = expression(k, j) * ps[k] / sumx[j];
    }
  }, ncores);
  return pMatrix;
}


//' @import Rcpp
//' @import Matrix
//' @title rcpp_pStrata_parallel
//' @name rcpp_pStrata_parallel
//' @description computes the partial
//' transcriptome evolutionary index (TEI) values combined into strata
//' @return sparseMatrix
//' @param expression ExpressionSet as sparseMatrix
//' @param ps named Phylostratum
//' @param psgroup ordered unique Phylostratum
//' @param ncores number of cores
//' @examples
//' ## load example PhyloExpressionSetExample
//'
//' data("PhyloExpressionSetExample", package="myTAI")
//'
//' ## convert into sparseMatrix - rownames GeneID
//'
//' spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]),
//'     "sparseMatrix")
//' rownames(spmat) <- PhyloExpressionSetExample$GeneID
//'
//' ## create named Phylostratum vector
//'
//' ps <- setNames(PhyloExpressionSetExample$Phylostratum,
//'     PhyloExpressionSetExample$GeneID)
//' psgroup <- sort(unique(ps))
//'
//' ## get pStrata
//' rcpp_pStrata_parallel(spmat, ps, psgroup)
//' @author Kristian K Ullrich
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_pStrata_parallel(const arma::sp_mat& expression,
                                          Rcpp::NumericVector ps,
                                          Rcpp::NumericVector psgroup,
                                          int ncores = 1){
  std::vector< std::string > psnames =  ps.attr("names");
  int n_col = expression.n_cols;
  int n_row = expression.n_rows;
  int n_psgroup = psgroup.length();
  Rcpp::NumericMatrix sMatrix(n_psgroup, n_col);
  Rcpp::NumericVector sumx(n_col);
  std::unordered_map<double, double> index_psgroup;
  for (size_t p = 0; p < n_psgroup; p++) {
    index_psgroup[psgroup[p]] = p;
  }
  RcppThread::ProgressBar bar(n_col, 1);
  RcppThread::parallelFor(0, n_col, [&] (int j) {
    for (size_t i = 0; i < n_row; i++) {
      sumx[j] += expression(i, j);
    }
    for (size_t k = 0; k < n_row; k++) {
      sMatrix(index_psgroup[ps[k]], j) += expression(k, j) * ps[k] / sumx[j];
    }
  }, ncores);
  return sMatrix;
}



