#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  
  if (a > lambda) {
    return a - lambda;
  } else if (a < -lambda) {
    return a + lambda;
  } else {
    return 0.0;
  }
  
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  int n = Xtilde.n_rows;
  double rss = accu(square(Ytilde - Xtilde * beta)) / (2.0 * n);
  double penalty = lambda * accu(abs(beta));
  return rss + penalty;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
  int p = Xtilde.n_cols;
  int n = Xtilde.n_rows;
  
  arma::colvec beta = beta_start;
  arma::colvec r = Ytilde - Xtilde * beta;
  
  double f_old = lasso_c(Xtilde, Ytilde, beta, lambda);
  
  while(true){
    for(int j = 0; j < p; j++){
      double beta_j_old = beta(j);
      double a_j = arma::dot(Xtilde.col(j), Xtilde.col(j)) / n;
      double c_j = arma::dot(Xtilde.col(j), r) / n + beta(j) * a_j;
      beta(j) = soft_c(c_j / a_j, lambda / a_j);
      r -= Xtilde.col(j) * (beta(j) - beta_j_old);
    }
    
    double f_new = lasso_c(Xtilde, Ytilde, beta, lambda);
    if(std::abs(f_new - f_old) < eps) break;
    f_old = f_new;
  }
  
  return beta;
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
  int p = Xtilde.n_cols;
  int L = lambda_seq.n_elem;
  
  arma::mat beta_mat(p, L, arma::fill::zeros);
  arma::colvec beta = arma::zeros(p); // warm start
  
  for(int l = 0; l < L; l++){
    beta = fitLASSOstandardized_c(Xtilde, Ytilde, lambda_seq(l), beta, eps);
    beta_mat.col(l) = beta;
}
  return beta_mat;
}