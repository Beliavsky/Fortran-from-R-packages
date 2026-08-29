# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]


//' Product of two matrices
//' 
//' @param A1 A matrix.
//' @param A2 A matrix.
//' @return Computes A1 * A2.
//' 
// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_product(Rcpp::NumericMatrix A1, Rcpp::NumericMatrix A2) {
  arma::mat AA1 = Rcpp::as<arma::mat>(A1);
  arma::mat AA2 = Rcpp::as<arma::mat>(A2);
  return(Rcpp::wrap(AA1 * AA2));
}


//' Inverse of a matrix
//' 
//' @param A A matrix.
//' @return Inverse of A. 
//' 
// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_inverse(Rcpp::NumericMatrix A) {
  arma::mat AA = Rcpp::as<arma::mat>(A);
  return(Rcpp::wrap(inv(AA)));
}
