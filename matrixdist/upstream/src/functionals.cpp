#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Univariate case

//' Laplace transform of a phase-type distribution
//'
//' Computes the Laplace transform at \code{r} of a phase-type distribution with
//'  parameters \code{alpha} and \code{S}.
//'
//' @param r Vector of real values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @return Laplace transform at \code{r}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector ph_laplace(Rcpp::NumericVector r, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector laplace(r.size());

  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;

  arma::mat aux_mat(1,1);

  arma::mat identity_matrix;
  identity_matrix.eye(size(S));

  for (int i{0}; i < r.size(); ++i) {
    aux_mat = alpha.t() * inv(identity_matrix * r[i] +  S * (-1.0)) * exit_vect;
    laplace[i] = aux_mat(0,0);
  }
  return laplace;
}

//' Pgf of a discrete phase-type distribution
//'
//' Computes the pgf at \code{z} of a discrete phase-type distribution with
//'  parameters \code{alpha} and \code{S}.
//'
//' @param z Vector of real values.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-transition matrix.
//' @return Laplace transform at \code{r}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector dph_pgf(Rcpp::NumericVector z, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector pgf(z.size());
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = e - (S * e);
  
  arma::mat aux_mat(1,1);
  
  arma::mat identity_matrix;
  identity_matrix.eye(size(S));
  
  for (int i{0}; i < z.size(); ++i) {
    aux_mat = alpha.t() * inv(identity_matrix / z[i] +  S * (-1.0)) * exit_vect;
    pgf[i] = aux_mat(0,0);
  }
  return pgf;
}

// Multivariate case

//' Bivariate phase-type joint Laplace
//'
//' @param r Matrix of values.
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @return Joint laplace at \code{r}.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector bivph_laplace(Rcpp::NumericMatrix r, arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22) {
  long n{r.nrow()};

  Rcpp::NumericVector laplace(n);

  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = (S22 * (-1)) * e;

  arma::mat identity_matrix1;
  identity_matrix1.eye(size(S11));

  arma::mat identity_matrix2;
  identity_matrix2.eye(size(S22));

  arma::mat aux_mat(1,1);

  for (int k{0}; k < n; ++k) {
    aux_mat = alpha.t() * inv(identity_matrix1 * r(k,0) +  S11 * (-1.0)) * S12 * inv(identity_matrix2 * r(k,1) +  S22 * (-1.0)) * exit_vect;
    laplace[k] = aux_mat(0,0);
  }
  return laplace;
}

