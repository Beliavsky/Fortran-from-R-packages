#ifndef transformations
#define transformations

#include <RcppArmadillo.h>

int n_pos (const arma::vec R);

arma::vec plus_states (const arma::vec R);

Rcpp::List tvr_ph(arma::vec alpha, arma::mat S, arma::vec R);

#endif
