//#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <math.h>
// // [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;

/* This whole 'permut' function has been adapted and taken from:
   http://gallery.rcpp.org/articles/stl-random-shuffle/

   - wrapper around R's RNG such that we get a uniform distribution over
   - [0,n) as required by the STL algorithm

*/
inline int randWrapper(const int& n)
{ 
  return floor(unif_rand() * n); 
}

NumericVector permut(const NumericVector& a)
{

    // already added by sourceCpp(), but needed standalone
    RNGScope scope;             

    // clone a into b to leave a alone
    NumericVector b = clone(a);

    random_shuffle(b.begin(), b.end(), randWrapper);

    return b;
}


// @export
// [[Rcpp::export]]
 NumericVector cpp_TAI(const NumericMatrix& ExpressionSet, const NumericVector& Phylostratum)
 {
 
    int nCols = ExpressionSet.ncol();
    int nRows = ExpressionSet.nrow();
    NumericVector results(nCols);
    for(int stage = 0; stage < nCols; stage++) {
	    double numerator = 0, divisor = 0;
	      for(int gene = 0; gene < nRows; gene++) {
		      numerator+= (double) Phylostratum[gene] * ExpressionSet(gene, stage);
		      divisor  += ExpressionSet(gene, stage);
	       }

	       results[stage] = numerator/divisor;
    }
    
    return results;

  }

// @export
// [[Rcpp::export]]
NumericMatrix cpp_bootMatrix(const NumericMatrix& ExpressionMatrix, const NumericVector& AgeVector, const int& permutations)
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
  NumericMatrix cpp_pMatrix(const NumericMatrix& ExpressionSet,const NumericVector& AgeVector)
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


