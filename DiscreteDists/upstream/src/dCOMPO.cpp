#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <cmath>
#include <complex>

// [[Rcpp::export]]
double z_cpp(double lambda, double nu, int max_terms = 1000,
             double tol = 1e-10) {
  double sum = 0.0;
  double term;
  int j = 0;

  while (j < max_terms) {
    double denom = std::pow(std::tgamma(j + 1), nu);  // (j!)^ν usando gamma(j+1)
    term = std::pow(lambda, j) / denom;
    sum += term;

    if (term < tol) break;  // termina si el termino es muy pequeno

    ++j;
  }

  return sum;
}

// [[Rcpp::export]]
NumericVector z_vec_cpp(NumericVector mu, NumericVector sigma) {

  int n = mu.size();
  NumericVector out(n);

  if (mu.size() != sigma.size()) {
    throw std::invalid_argument("Error: Vectors mu and sigma must be of the same length.\n");
    return out;
  }

  for(int i = 0; i < n; ++i) {
    out[i] = z_cpp(mu[i], sigma[i]);
  }
  return out;
}

// [[Rcpp::export]]
 double d1_dldm_compo_cpp(double lambda, double nu,
                          int max_terms = 1000, double tol = 1e-10) {
   double sum = 0.0;
   double term;
   int j = 1;

   while (j < max_terms) {
     double numer = j * std::pow(lambda, j-1);
     double denom = std::pow(std::tgamma(j + 1), nu);
     term = numer / denom;
     sum += term;

     if (term < tol) break;  // termina si el termino es muy pequeno

     ++j;
   }

   return sum;
 }


// [[Rcpp::export]]
NumericVector d1_vec_dldm_compo_cpp(NumericVector mu, NumericVector sigma) {

  int n = mu.size();
  NumericVector out(n);

  if (mu.size() != sigma.size()) {
    throw std::invalid_argument("Error: Vectors mu and sigma must be of the same length.\n");
    return out;
  }

  for(int i = 0; i < n; ++i) {
    out[i] = d1_dldm_compo_cpp(mu[i], sigma[i]);
  }
  return out;
}

// [[Rcpp::export]]
 double d2_dldd_compo_cpp(double lambda, double nu,
                          int max_terms = 1000, double tol = 1e-10) {
   double sum = 0.0;
   double term;
   int j = 2;

   while (j < max_terms) {
     double fact = std::tgamma(j + 1); // j! = tgamma(j+1)
     double numer = std::pow(lambda, j) * std::log(fact);
     double denom = std::pow(fact, nu);
     term = numer / denom;
     sum += term;

     if (term < tol) break;  // termina si el termino es muy pequeno

     ++j;
   }

   return sum;
}

// [[Rcpp::export]]
NumericVector d2_vec_dldd_compo_cpp(NumericVector mu, NumericVector sigma) {

  int n = mu.size();
  NumericVector out(n);

  if (mu.size() != sigma.size()) {
    throw std::invalid_argument("Error: Vectors mu and sigma must be of the same length.\n");
    return out;
  }

  for(int i = 0; i < n; ++i) {
    out[i] = d2_dldd_compo_cpp(mu[i], sigma[i]);
  }
  return out;
}

// [[Rcpp::export]]
 double dCOMPO_single(double x, double mu=1, double sigma=1, bool log=false) {
   if (sigma <= 0 || mu <= 0) {
     throw std::runtime_error("parameter sigma and mu must be positive!");
   }

   double res;

   if (x < 0) {
     res = std::log(0);
   }
   else {
     double p1 = x * std::log(mu) - sigma * std::lgamma(x + 1);
     double p2 = z_cpp(mu, sigma);
     res = p1 - std::log(p2);
   }

   if (log) {
     return res;
   } else {
     return std::exp(res);
   }
 }

// [[Rcpp::export]]
 NumericVector dCOMPO_vec(NumericVector x, NumericVector mu,
                            NumericVector sigma, LogicalVector log) {
   int n = x.size();
   NumericVector out(n);

   for(int i = 0; i < n; ++i) {
     out[i] = dCOMPO_single(x[i], mu[i], sigma[i], log[i]);
   }
   return out;
 }

