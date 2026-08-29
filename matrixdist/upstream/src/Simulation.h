#include <RcppArmadillo.h>

arma::mat embedded_mc(arma::mat S);

arma::mat cumulate_matrix(arma::mat A);

arma::vec cumulate_vector(arma::vec A);

long initial_state(arma::vec cum_alpha, double u);

long new_state(long prev_state, arma::mat cum_embedded_mc, double u);

Rcpp::NumericVector rphasetype(int n, arma::vec alpha, arma::mat S);

Rcpp::NumericVector riph(int n, Rcpp::String dist_type, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta);

Rcpp::NumericVector rmatrixgev(int n, arma::vec alpha, arma::mat S, double mu, double sigma, double xi = 0);
