// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include "bzinb.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
using namespace std;
using namespace Rcpp;
//digamma overflows when arg >> 600. So approximate it with log(x - 0.5). See https://en.wikipedia.org/wiki/Digamma_function#Computation_and_approximation
// digamma(x) is in (log(x - 0.5), log(x)) and is much closer to the left limit.
#define digamma2(a) (((a) < (600))?(boost::math::digamma(x)):(log(x - 0.5)))

#define EPSILON1 1e-9  //for inverse-gamma
#define EPSILON2 1e-9  //for M-step
// #define DEBUG
// #define DEBUG2
// #define DEBUG3
//bool debug = true;

// 2. optimization
double inv_digamma2(double x, double y) 
{ 
  double h = (digamma2(x) - y) / boost::math::trigamma(x);
  while (fabs(h) >= EPSILON1)
  {
    h = (digamma2(x) - y) / boost::math::trigamma(x);
    while (h > x) {
      #ifdef DEBUG6 
        Rcout <<"inv_digamm hit zero:" << h << " ";
      #endif
      h = x/2.0;  // use x/2.0 (bisection)
    }
// #ifdef DEBUG 
//   Rcout << "x:" << x << ", y:" << y << ", h:" << h << ", digamma(x): " << digamma2(x) <<  ", tri(x):" << boost::math::trigamma(x) << " " << endl;
// #endif  
    x = x - h;
  }
  // if (debug) {Rcout << endl << "inv_digamm = " << x << endl;}
  #ifdef DEBUG 
    Rcout << "inv_digamm = " << x << " ";
  #endif
  
  return(x);
}

// [[Rcpp::export]]
double inv_digamma(double x, double y) 
{ 
  if (x < 600) {
    return(inv_digamma2(x, y));
  } else {
    return(exp(x + 0.5));
  }
}

// optimization


void inv_digamma_vec(double lb[1], NumericVector &expt, Rcpp::NumericVector &a, double idgam[3]) 
{ 
  //double idgam[3];
  //double *idgam = new double[3];
#ifdef DEBUG
  Rcout << endl<< "idgam[0] start. expt[4] = " << expt[4] << ", lb[0] = " << lb[0] << endl;
#endif
  idgam[0] = inv_digamma(a[0], expt[4] - lb[0]);
#ifdef DEBUG
  Rcout  << endl<< "idgam[0] = " << idgam[0] << endl <<"idgam[1] start. expt[5] = " << expt[5] << ", lb[0] = " << lb[0] << endl;
#endif
  idgam[1] = inv_digamma(a[1], expt[5] - lb[0]);
#ifdef DEBUG
  Rcout  << endl<< "idgam[1] = " << idgam[1] << endl <<"idgam[2] start. expt[6] = " << expt[6] << ", lb[0] = " << lb[0]  << endl;
#endif
  idgam[2] = inv_digamma(a[2], expt[6] - lb[0]);
#ifdef DEBUG
  Rcout  << endl<< "idgam[2] = " << idgam[2] << endl;
#endif
  //return(idgam);
}

// double objective(double lb, NumericVector &expt, NumericVector &a, NumericVector &idgam) 
// { 
//   inv_digamma_vec(lb, expt, a, idgam);
//   // Rcout << "idgam123: "<< idgam[0] << " "<< idgam[1] << " "<< idgam[2] << endl;
//   //  inv_digamma_vec(b1, expt, a, idgam);
//   // double result = (double);
//   return (log(idgam[0] + idgam[1] + idgam[2]) + lb - log(expt[1] + expt[2] + expt[3]));
// }
// 
// // Derivative
// double derivFunc(double lb, NumericVector &expt, NumericVector &a, NumericVector &idgam)
// {
//   inv_digamma_vec(lb, expt, a, idgam);
//   double result = 0.0;
//   for (int i = 0; i < 3; i++) {
//     result += 1/ boost::math::trigamma(idgam[i]);
//   }
//   return ((result)/(idgam[0] + idgam[1] + idgam[2]) + 1.0);
// }

// Derivative
double hfunc(double lb[1], NumericVector &expt, NumericVector &a, double idgam[3])
{
  double result = 0.0;
  double obj;
  
  if (((expt[4] - lb[0] < 600) && (expt[5] - lb[0] < 600)) && (expt[6] - lb[0] < 600)) {
    // when everything is not large.
    inv_digamma_vec(lb, expt, a, idgam);
    obj = (log(idgam[0] + idgam[1] + idgam[2]) + lb[0] - log(expt[1] + expt[2] + expt[3]));
    
  #ifdef DEBUG 
    Rcout << "objective = " << obj << " ";
  #endif
    for (int i = 0; i < 3; i++) {
      // Rcout << "idgam[0:2]" << idgam[0] << " " << idgam[1] << " " << idgam[2] << endl;
      result += (1/ boost::math::trigamma(idgam[i]));
    }
    result = - result / (idgam[0] + idgam[1] + idgam[2]) + 1.0;
    result = obj / result;
    return (result);
  } else {
    // approximation for large values.
    // idgam will not be updated.
    Rcout << "large idgamma activated." << endl;
    double max_expt = std::max(std::max(expt[4], expt[5]), expt[6]) - 0.5;
    obj = max_expt + log(exp(expt[4] - max_expt) + exp(expt[5] - max_expt) + exp(expt[6] - max_expt)) -
      log(expt[1] + expt[2] + expt[3]);
    
  #ifdef DEBUG 
    Rcout << "objective = " << obj << " ";
    Rcout << "idgam[0:2]" << idgam[0] << " " << idgam[1] << " " << idgam[2] << endl;
    Rcout << "expt[4:6] - lb = " << expt[4] - lb[0] << " " <<  expt[5] - lb[0] << " " <<  expt[6] - lb[0] << endl;
  #endif
    result = obj * ( 1 + 2 * exp(max_expt + 0.5) - lb[0] + 0.5);
    return (result);
  }
  

}

void opt_lb(double lb[1], NumericVector &expt, NumericVector &a, double idgam[3])
{
  int i = 1;
  //double* lb = log(b1);
  // double h =  objective(lb, expt, a, idgam) / derivFunc(lb, expt, a, idgam);
  double h = hfunc(lb, expt, a, idgam);
  // if (debug) {Rcout << "h = " << h << " ";}
#ifdef DEBUG 
  Rcout << "h = " << h << " ";
#endif
#ifdef DEBUG3
  Rcout << "opt_iter = 1 ";  
#endif
  while ((fabs(h) >= EPSILON2) & (i < 11))  //sometimes h oscilliates around zero 
                // (usually close to solution with only a few iterations, as the initial value is close to truth.)
                // See p.10 of http://mathfaculty.fullerton.edu/mathews//n2003/newtonsmethod/Newton'sMethodProof.pdf
  {
    i += 1;
    // h = objective(lb, expt, a, idgam)/derivFunc(lb, expt, a, idgam);
    // while (h > b1) {
    //   // Rcout << "opt_b1 hit zero. ";
    //   h /= 2.0;
    // }
    // Rcout << "h = " << h << " ";
#ifdef DEBUG3 
    Rcout << "(i= " << i<< ", lb= " << lb[0] << ", h= " << h << ", lb-h=" << lb[0]-h << ").. ";
#endif
    lb[0] -= h;
    h = hfunc(lb, expt, a, idgam);
  }
#ifdef DEBUG2 
  Rcout << "lb =" << lb[0]  << ", b1 = " << exp(lb[0]) << " / ";
#endif
  
// if (debug) {Rcout << "lb =" << lb  << ", b1 = " << exp(lb) << endl;}
// if (debug) {Rcout << "lb =" << lb  << endl;}
  // double* idgam = inv_digamma_vec(lb, expt, a);
  a[0] = idgam[0];
  a[1] = idgam[1];
  a[2] = idgam[2];
  // Rcout << "a123: "<< a[0] << " "<< a[1] << " "<< a[2] << " " << exp(lb) << endl;
}
