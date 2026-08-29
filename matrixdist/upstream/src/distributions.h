#include <RcppArmadillo.h>

Rcpp::NumericVector phdensity(Rcpp::NumericVector x, arma::vec alpha, arma::mat S);

Rcpp::NumericVector phcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, bool lower_tail = true);

Rcpp::NumericVector mweibullden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta);

Rcpp::NumericVector mweibullcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta, bool lower_tail = true);

Rcpp::NumericVector mparetoden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta) ;

Rcpp::NumericVector mparetocdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta, bool lower_tail = true);
  
Rcpp::NumericVector mlognormalden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta);

Rcpp::NumericVector mlognormalcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta, bool lower_tail = true);

Rcpp::NumericVector mloglogisticden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta);

Rcpp::NumericVector mloglogisticden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta);

Rcpp::NumericVector mgompertzden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta);

Rcpp::NumericVector mgompertzcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, double beta, bool lower_tail = true);

Rcpp::NumericVector mgevden(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta);

Rcpp::NumericVector mgevcdf(Rcpp::NumericVector x, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta, bool lower_tail = true);
