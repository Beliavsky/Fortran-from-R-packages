#include <Rcpp.h>
using namespace Rcpp;


// Function rowSumsModN_cpp:
// Takes the vector x and gives the same output as:
// rowSums(matrix(x, ncol=N, nrow=length(x)/N, byrow =TRUE)) in R.
// Used to compute the time-t posterior density by summing over
// a single x_t value (and across all other latent factors in the DNF function).
// Inputs:
//  - x : Vector of probabilities.
// 	- N : Length of rows for the matrix that we want to sum over.
// Output:
// 	- ans : Result of rowSums(matrix(x, ncol=N, nrow=length(x)/N, byrow =TRUE))

// [[Rcpp::export]]
NumericVector Cpp_rowSums_modN(NumericVector x, int N) {
  int sizeX = x.size();
  NumericVector ans(N);
  for (int j = 0; j < N; j++) {
    for (int i = j; i < sizeX; i = i + N) {
      ans[j] += x[i];
    }
  }
  return ans;
}

// Function prodFun_cpp:
// This function multiplies the elements of the filtering distribution with there corresponding
// x_tmin1 entries from the expand.grid in the probs vector
// (the products in Equation 8 of Bégin and Boudreault 2020).
// Inputs:
//  - probs : Product of transition and measurement probabilities.
// 	- filter : x_tmin1 filtering distribution estimate.
// Output:
// 	- probs : Product of transition, measurement probabilities and the appropriate
// x_tmin1 filtering distribution estimates.

// [[Rcpp::export]]
NumericVector Cpp_prodfun(NumericVector probs, NumericVector filter) {
  int n1 = filter.size(); // length of the filtering grid, N
  int n2 = probs.size(); // length of the expand.grid for all latent factor, N*N*K*R
  int n3 = n2/n1; // N*K*R
  double x0;
  for (int j = 0; j < n3; j++) {
    x0 = filter[j % n1]; //Each filter elements are repeated N times in every N*N vector elements
    for (int i = j*n1; i < (j+1)*n1; i++) {
     probs[i] *= x0;
    }
  }
  return probs;
}


// Function dnormProd_cpp:
// This function computes the normal density associated with
// the measurement equation (r() in Bégin and Boudreault 2020))
// and multiplies them, element by element, with the elements of the vector p1.
// Here, p1 will be the product of the other factors' densities.
// Inputs:
//  - rets : Product of transition and measurement probabilities.
// 	- mu : x_tmin1 filtering distribution estimate.
// 	- sigma : x_tmin1 filtering distribution estimate.
// 	- p : x_tmin1 filtering distribution estimate.

// Output:
// 	- probs : Product of transition, measurement probabilities and the appropriate
// x_tmin1 filtering distribution estimates.

// [[Rcpp::export]]
NumericMatrix dnorm_cpp_prod(NumericVector rets, NumericVector mu, NumericVector sig, NumericVector p){
  double c = 1/sqrt(2*3.14159265358979);
  int NNKR = mu.size(); // Size of the expand.grid for all our factors (N*N*K*R)
  int t = rets.size(); // x is our vector of returns and x.size = t
  NumericMatrix ans(NNKR, t);
  double x0, s0, x1;
  for (int j = 0; j < t; j++) {
    x1 = rets[j]; // Returns observation at each time
    for (int i = 0; i < NNKR; i++) {
      // Computing the normal densities for each combination in the NNKR grid at time j.
      x0 = x1-mu[i];
      s0 = sig[i];
      ans(i,j) =(exp(-x0*x0/(2*s0*s0))*c/s0)*p[i];
    }
  }
  return ans;
}

