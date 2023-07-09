//#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
#include <math.h>
#include <map>
#include <random>
#include <RcppThread.h>
#include <string.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
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

// Initilizing the random number generator outside of the function
std::random_device rng;
std::mt19937_64 urng(rng());
std::linear_congruential_engine<unsigned int, 48271, 1, 65536> lcg(rng());

Eigen::MatrixXd permut_mat(const Eigen::VectorXd& a,const int& permutations) {
  // already added by sourceCpp(), but needed standalone
  Eigen::MatrixXd permutedMat(permutations, a.size());
  Eigen::VectorXd shuffledVec = a;
  
  
  for (int i = 0; i < permutations; i++) {
    std::shuffle(shuffledVec.data(), shuffledVec.data() + shuffledVec.size(), lcg);
    permutedMat.row(i) = shuffledVec.transpose();
  }
  return permutedMat;
}

NumericVector permut(const NumericVector& a)
{
  
  // already added by sourceCpp(), but needed standalone
  // RNGScope scope;             
  
  // clone a into b to leave a alone
  NumericVector b = clone(a);
  
  std::shuffle(b.begin(), b.end(), lcg);
  
  return b;
}

// @export
// [[Rcpp::export]]
NumericMatrix cpp_bootMatrix_old(const NumericMatrix& ExpressionMatrix, const NumericVector& AgeVector, const int& permutations)
{
  
  int nCols = ExpressionMatrix.ncol();
  int nRows = ExpressionMatrix.nrow();
  NumericVector Divisor(nCols);
  NumericMatrix fMatrix(nRows,nCols);
  NumericVector sampledVector(AgeVector.size());
  NumericVector AgeValues(nCols);
  NumericMatrix bootM(permutations,nCols);
  
  for(int j = 0; j < nCols; j++){
    double div = 0;
    for(int i = 0; i < nRows; i++){
      div += ExpressionMatrix(i,j);    
    }
    
    Divisor[j] = div;
    
  }
  
  for(int l = 0; l < nCols; l++){
    for(int k = 0; k < nRows; k++){
      fMatrix(k,l) = ExpressionMatrix(k,l)/Divisor[l];    
    }
  }
  
  
  for(int n = 0; n < permutations; n++){
    
    //sampledVector = Rcpp::RcppArmadillo::sample(AgeVector, AgeVector.size(), FALSE);
    sampledVector = permut(AgeVector);
    
    for(int b = 0; b < nCols; b++){
      double total = 0;
      for(int a = 0; a < nRows; a++){
        total += (sampledVector[a] * fMatrix(a,b));    
      }
      
      bootM(n,b) = total;
      
    }
  }
  
  return bootM;
}





// @export
// [[Rcpp::export]]
Eigen::VectorXd cpp_TAI(const Eigen::MatrixXd& ExpressionMatrix, const Eigen::VectorXd& Phylostratum) {
  
  Eigen::VectorXd Divisor = ExpressionMatrix.colwise().sum();
  Eigen::MatrixXd fMatrix = ExpressionMatrix.array().rowwise() / Divisor.transpose().array();
  Eigen::VectorXd total = Phylostratum.transpose() * fMatrix;
  
  return total;
}
// @export
// [[Rcpp::export]]
Eigen::MatrixXd cpp_bootMatrix(const Eigen::MatrixXd& ExpressionMatrix, const Eigen::VectorXd& AgeVector, const int& permutations) 
{
        Eigen::VectorXd Divisor = ExpressionMatrix.colwise().sum();
        Eigen::MatrixXd fMatrix = ExpressionMatrix.array().rowwise() / Divisor.transpose().array();
        Eigen::MatrixXd permMatrix = permut_mat(AgeVector,permutations);
        Eigen::MatrixXd bootM = permMatrix * fMatrix;
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
NumericMatrix cpp_pMatrix_old(const NumericMatrix& ExpressionSet,const NumericVector& AgeVector)
{
  
  int nRows = ExpressionSet.nrow();
  int nCols = ExpressionSet.ncol(); 
  NumericMatrix results(nRows,nCols);
  NumericVector DivisorVector(nCols);
  
  for(int stage = 0; stage < nCols; stage++) {
    double divisor = 0;
    for(int gene = 0; gene < nRows; gene++) {
      divisor  += ExpressionSet(gene, stage);
    }
    
    DivisorVector[stage] = divisor;
  }
  
  for(int stage = 0; stage < nCols; stage++){
    for(int gene = 0; gene < nRows; gene++){
      results(gene,stage) = (double) AgeVector[gene] * (ExpressionSet(gene, stage)/DivisorVector[stage]);
    }
  }
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

// @export
// [[Rcpp::export]]
NumericMatrix cpp_omitMatrix(const NumericMatrix& ExpressionSet, const NumericVector& AgeVector){
  
  int nRows = ExpressionSet.nrow();
  int nCols = ExpressionSet.ncol(); 
  NumericMatrix ResultMatrix(nRows,nCols);
  NumericVector NumeratorVector(nCols);
  NumericVector DivisorVector(nCols);
  
  
  for(int stage = 0; stage < nCols; stage++) {
    double numerator = 0, divisor = 0;
    for(int gene = 0; gene < nRows; gene++) {
      numerator+= (double) AgeVector[gene] * ExpressionSet(gene, stage);
      divisor  += ExpressionSet(gene, stage);
    }
    
    NumeratorVector[stage] = numerator;
    DivisorVector[stage] = divisor;
  }
  
  for(int stage = 0; stage < nCols; stage++){
    double newNumerator = 0, newDivisor = 0;
    for(int gene = 0; gene < nRows; gene++){
      newNumerator = (double) NumeratorVector[stage] - (AgeVector[gene] * ExpressionSet(gene, stage));
      newDivisor = (double) DivisorVector[stage] - ExpressionSet(gene, stage);
      ResultMatrix(gene,stage) = newNumerator / newDivisor;
    }
  }
  
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



