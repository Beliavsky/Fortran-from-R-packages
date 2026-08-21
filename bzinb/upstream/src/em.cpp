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
#define ITER_ALLOWANCE 100   // number of iterations allowed after finding the peak
// #define DEBUG4
// #define DEBUG5
// #define DEBUG7

// 3. EM
// [[Rcpp::export]]
List em(NumericVector& param2, IntegerVector &xvec, IntegerVector &yvec, 
        IntegerVector &freq, int &n, 
        int &se, int &maxiter, double &tol, int showFlag,
        int bnb)
{
  NumericVector param = clone(param2);
  double param_diff = 1.0;
  double param_old[9];
  double idgam[3];
  double lb[1];
  NumericVector expt(12, 0.0);
  NumericVector s_i(8, 0.0);
  NumericVector info(1 + se * 63L, 0.0);  // if se = 1, 8 x 8 matrix, o/w a scalar zero.
  IntegerVector iter(1, 1L);
  IntegerVector nonconv(1, 0L);
  NumericVector trajectory(maxiter + 1L, 0.0); 
  
  // storage for max_likelihood.
  double param_max[9];
  int iter_max = 0;
  double expt_max[12];
  
  for (int i = 0; i < 9; i++) {
    param[i] = (double) param2[i];  //initializing param2 with param
  }
  
  //cout << "maxiter = " << maxiter << " iter = " << iter[0] << endl;  
  while(maxiter >= iter[0] && param_diff > tol)
  {
    if (iter[0] < 3) {
      expt_max[0] = expt[0];
    }
    //cout << "maxiter = " << maxiter << " iter1 = " << iter[0] << endl;  
    for(int i = 0;i < 9;i++)
    {
      param_old[i] = param[i];
      // cout << "param [" << i << "] = " << param[i] << " ";
    }
    dBvZINB_Expt_vec(xvec, yvec, freq, n, param[0], param[1], param[2], param[3], param[4], 
                     param[5], param[6], param[7], param[8], expt, s_i, info, 0, bnb);
    
    if ((showFlag > 0) & (iter[0] >= showFlag)) {
      Rcout << "iter " << iter[0]  << " lik " << expt[0] <<
        ", a0 " << param[0] << ", a1 " << param[1] << ", a2 " << param[2] <<
          ", b1 " << param[3] << ", b2 " << param[4];
      if (!bnb) {
        Rcout << ", p1 " << param[5] << ", p2 " << param[6] << ", p3 " << param[7] << 
          ", p4 " << param[8];
      }
      Rcout << endl;
    }
    
#ifdef DEBUG4
    for (int i = 0; i < 8; i++) 
    {
      if (i == 0) {Rcout << "param: ";}
      Rcout << param[i] << " ";
      if (i == 7) {Rcout <<  endl;}
    }
    for (int i = 0; i < 12; i++) 
    {
      if (i == 0) {Rcout << "expt: ";}
      Rcout << expt[i] << " ";
    }
    Rcout << endl;
#endif
    
    for(int i = 0;i < 4; i++)
    {
      param[5 + i] = expt[7 + i];
    }
    lb[0] = log(param[3]);
    
#ifdef DEBUG5
    for (int i = 0; i < 3; i++) 
    {
      if (i == 0) {Rcout << "idgam: ";}
      Rcout << idgam[i] << " ";
    }
    Rcout <<  "lb = " << lb[0]  << endl;
    Rcout << "before opt_lb! (of iter" << iter << ") ";
#endif
    //if (iter[0] > 2058) {cout << "before opt_lb" << endl;}
    
    // Finding optimized a0, a1, a2, b1
#ifdef DEBUG7
    Rcout << "starting opt_lb. lb = " << lb[0] << endl;
#endif  
    // after finding "nice" initial lb, run opt_lb.
    opt_lb(lb, expt, param, idgam);
#ifdef DEBUG5
  Rcout << ", after opt_lb! (of iter" << iter << ") "<< lb[0] << endl;
#endif
    //if (iter[0] > 2058) {cout << "before opt_lb" << endl;}
    double delta = expt[11]*1.0 / (expt[1] + expt[3]);
    param[3] = exp(lb[0]);
    param[4] = param[3] * delta;
    
    param_diff = 0.0;
    
    for(int i = 0;i < 9;i++)
    {
      double dif = fabs(0.0 + param[i] - param_old[i]);
      if( dif > param_diff) param_diff = dif;
    }
    
    // updating max lik
    if (expt[0] > expt_max[0]) {
      iter_max = iter[0];
      for(int i = 0;i < 9;i++) {
        param_max[i] = param_old[i];
      }
      for(int i = 0;i < 12;i++) {
        expt_max[i] = expt[i];
      }
    }
    
    //If iteration lasts more than ITER_ALLOWANCE times without improvement of lik, stops.
    if (iter[0] >= iter_max + ITER_ALLOWANCE) { 
      break;
    }
    trajectory[iter[0]] = expt[0];
    iter[0] += 1;
  }
  
  if (maxiter <= iter[0] &&  param_diff > tol) {nonconv[0] = 1;}  
  
  // cout << "pluging in max_lik iter = " << iter_max << " lik = " << expt_max[0] << " instead of lik = " << expt[0] << endl;
  // replacing back final result with historical maximum
  iter[0] = iter_max;
  for(int i = 0;i < 9;i++) {
    param[i] = param_max[i];
  }
  for(int i = 0;i < 12;i++) {
    expt[i] = expt_max[i];
  }

  if (se == 1) { //updating expt and calculate SE when called for.
    dBvZINB_Expt_vec(xvec, yvec, freq, n, param[0], param[1], param[2], param[3], param[4], 
                     param[5], param[6], param[7], param[8], expt, s_i, info, se, bnb);
  }

  // cout << "a " << param[0] << " " << param[1] << " " << param[2] << " b " << param[3] << " " << param[4] << " pi "
  //      << param[5] << " " << param[6] << " "  << param[7] << " " << param[8] << endl;
  // for (int i = 0; i < 9; i++) {
  //   param2[i] = param[i];  //returning param to param2
  // }
  
  // List z = List::create(param, xvec, yvec, freq, n, expt, info, se, 
  //                       iter, nonconv, trajectory, bnb);
  List z = List::create(param, expt, info, iter, nonconv, trajectory);
  
  return z;
}
