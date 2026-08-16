#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <cmath>
#include <complex>

// [[Rcpp::export]]
 double f11_cpp(double gamma, double lambda,
                int maxiter_series = 10000,
                double tol = 1.0e-10) {
   double fac  = 1.0;
   double temp = 1.0;
   double L    = gamma;
   double series = temp;
   double f11;

   for (int n = 1; n <= maxiter_series; ++n) {
     fac = fac * lambda / L;
     series = temp + fac;

     if (std::abs(series - temp) < tol) {
       f11 = series;  // Assuming series is already real in this context
       return f11;
     }

     temp = series;
     L += 1;
   }

   f11 = series;  // Assuming series is already real in this context
   return f11;
 }

// [[Rcpp::export]]
 double dHYPERPO_single(double x, double mu=1, double sigma=1, bool log=false) {
   if (sigma <= 0 || mu <= 0) {
     throw std::runtime_error("parameter sigma and mu must be positive!");
   }

   double res;

   if (x < 0) {
     res = std::log(0);
   }
   else {
     double p1 = x * std::log(mu) - lgamma(sigma + x) + lgamma(sigma);
     double temp_f11 = f11_cpp(sigma, mu);
     double p2 = std::log(temp_f11);
     res = p1 - p2;
   }

   if (log) {
     return res;
   } else {
     return std::exp(res);
   }
 }

// [[Rcpp::export]]
 NumericVector dHYPERPO_vec(NumericVector x, NumericVector mu,
                            NumericVector sigma, LogicalVector log) {
   int n = x.size();
   NumericVector out(n);

   for(int i = 0; i < n; ++i) {
     out[i] = dHYPERPO_single(x[i], mu[i], sigma[i], log[i]);
   }
   return out;
 }


// [[Rcpp::export]]
double fun_exp_6_single_cpp(double x, double media, double gamma) {
  double res;
  res = x - (gamma-1) * (1-1/f11_cpp(gamma, x)) - media;
  return res;
}

// [[Rcpp::export]]
NumericVector fun_exp_6_vec_cpp(NumericVector x,
                                     NumericVector media,
                                     NumericVector gamma) {

  int n = media.size();
  NumericVector out(n);

  for(int i = 0; i < n; ++i) {
    out[i] = fun_exp_6_single_cpp(x[i], media[i], gamma[i]);
  }
  return out;
}


// [[Rcpp::export]]
double obtaining_lambda_single_cpp(double media,
                                   double gamma,
                                   double tol = 1e-10, int max_iter = 1000) {

  double lower = std::min(media, std::max(media + gamma - 1, gamma * media));
  double upper = std::max(media, std::min(media + gamma - 1, gamma * media));

  double f_lower = fun_exp_6_single_cpp(lower, media, gamma);
  double f_upper = fun_exp_6_single_cpp(upper, media, gamma);

  if (f_lower * f_upper > 0.0) {
    throw std::invalid_argument("Function has the same sign at the endpoints.");
  }

  double mid;
  double f_mid;

  if (gamma == 1) {
    mid = media;
    return mid;
  }

  else {
    for (int i = 0; i < max_iter; ++i) {
      mid = 0.5 * (lower + upper);
      f_mid = fun_exp_6_single_cpp(mid, media, gamma);

      if (std::fabs(f_mid) < tol || (upper - lower) / 2 < tol) {
        return mid;
      }

      if (f_lower * f_mid < 0) {
        upper = mid;
        f_upper = f_mid;
      } else {
        lower = mid;
        f_lower = f_mid;
      }
    }
  }

  return mid;

  throw std::runtime_error("Maximum iterations exceeded without convergence.");
}

// [[Rcpp::export]]
NumericVector obtaining_lambda_vec_cpp(NumericVector media,
                                       NumericVector gamma) {

  int n = media.size();
  NumericVector out(n);

  for(int i = 0; i < n; ++i) {
    out[i] = obtaining_lambda_single_cpp(media[i], gamma[i]);
  }
  return out;
}
