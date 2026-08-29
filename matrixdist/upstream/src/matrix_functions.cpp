#include <Rcpp.h>
using namespace Rcpp;


//' Default size of the steps in the RK
//' 
//' Computes the default step length for a matrix \code{S} to be employed in the
//'  RK method.
//'  
//' @param S Sub-intensity matrix.
//' @return The step length for \code{S}.
//' 
// [[Rcpp::export]]
double default_step_length(const NumericMatrix & S) {
  double h{-0.1 / S(0,0)};
  
  for (int i{1}; i < S.nrow(); ++i) {
    if (h > -0.1 / S(i,i)) {
      h = -0.1 / S(i,i);
    }
  }
  return h;
}


//' Applies the inverse of the GEV transformation but giving back the resulting 
//'  vector in reverse order
//' 
//' Used for EM step in RK.
//' 
//' @param obs The observations.
//' @param weights Weights of the observations.
//' @param beta Parameters of the GEV.
//' 
// [[Rcpp::export]]
List revers_data_trans(const NumericVector & obs, const NumericVector & weights, const NumericVector & beta) {
  int n{static_cast<int>(obs.size())};
  NumericVector trans_obs(n);
  NumericVector trans_weights(n);
  if (beta[2] == 0) { // Gumbel
    for (int i{0}; i < n; ++i) {
      trans_obs[i] = exp(-(obs[n - i - 1] - beta[0]) / beta[1]) ;
      trans_weights[i] = weights[n - i - 1];
    }
  }
  else { // GEVD
    for (int i{0}; i < n; ++i) {
      trans_obs[i] = pow(1 + beta[2] * (obs[n - i - 1] - beta[0]) / beta[1], -1 / beta[2]);
      trans_weights[i] = weights[n - i - 1];
    }
  }
  
  List L = List::create(Named("obs") = trans_obs, _["weight"] = trans_weights);
  return L;
}


//' Clone a vector 
//' 
//' @param v A vector.
//' @return A clone of the vector.
//' 
// [[Rcpp::export]]
NumericVector clone_vector(NumericVector v) {
  NumericVector new_v = clone(v);
  return new_v;
}


//' Clone a matrix 
//' 
//' @param m A matrix.
//' @return A clone of the matrix.
//' 
// [[Rcpp::export]]
NumericMatrix clone_matrix(NumericMatrix m) {
  NumericMatrix new_m = clone(m);
  return new_m;
}
