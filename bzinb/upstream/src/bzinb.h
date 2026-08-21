#ifndef BZINB_H
#define BZINB_H

#include <Rcpp.h>
#include <math.h>
#include <string>
#include <iostream>
#include <numeric>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
// 
// double EPSILON1;
// double EPSILON2;

void l1(int& x, int& y, double& a0, double& a1, double& a2, int &k, int& m, 
        double& result, double adjj = 0);
void l1_c (double& t1, double& t2, int& k, int& m, double& result, double adjj);
void l1_AC (double& t1, double& t2, int& x, int& y, double& a0, double& a1, 
            double& a2, int& k, int& m, double& result, double adjj = 0);
void l2_A (int& x, double& a0, double& a1, double& a2, int& k, 
           double& result, double adjj);
void l3_A (int& y, double& a0, double& a1, double& a2, int& m, 
           double& result, double adjj);
void R0_E1(int& x, int& y, int& k, int& m, double& a0, double& result);
double log_R0_E1(int& x, int& y, int& k, int& m, double& a0);
double	log_R0_E2(int& x, double& a0, int& k);
double	log_R0_E3(int& y, double& a0, int& m);
void R1_E1(int& k, double& a1, double& result);
double log_R1_E1(int& k, double& a1);

void dBvZINB_Expt(int &x, int &y, int &freq, double &a0, double &a1, double &a2,
                  double &b1, double &b2, double &p1, double &p2, 
                  double &p3, double &p4,
                  Rcpp::NumericVector &expt, Rcpp::NumericVector &s_i, Rcpp::NumericVector &info, 
                  int se, int bnb);
void dBvZINB_Expt_vec(Rcpp::IntegerVector &xvec, Rcpp::IntegerVector &yvec, 
                      Rcpp::IntegerVector &freq, 
                      int &n, double &a0, double &a1, double &a2,
                      double &b1, double &b2, double &p1, double &p2, 
                      double &p3, double &p4,
                      Rcpp::NumericVector &expt, Rcpp::NumericVector &s_i, Rcpp::NumericVector &info,
                      int se, int bnb);

void dBvZINB_Expt_direct(int &x, int &y, int &freq, double &a0, double &a1, double &a2,
                         double &b1, double &b2, double &p1, double &p2, 
                         double &p3, double &p4,
                         Rcpp::NumericVector &l_sum, Rcpp::NumericVector &s_i, 
                         Rcpp::NumericVector &info);
void dBvZINB_Expt_direct_vec(Rcpp::IntegerVector &xvec, Rcpp::IntegerVector &yvec, 
                             Rcpp::IntegerVector &freq, 
                             int &n, double &a0, double &a1, double &a2,
                             double &b1, double &b2, double &p1, double &p2, 
                             double &p3, double &p4,
                             Rcpp::NumericVector &l_sum, Rcpp::NumericVector &s_i, 
                             Rcpp::NumericVector &info);

double inv_digamma(double x, double y);
void inv_digamma_vec(double lb[1], Rcpp::NumericVector &expt, Rcpp::NumericVector &a, 
                     double idgam[3]);
double hfunc(double lb[1], Rcpp::NumericVector &expt, Rcpp::NumericVector &a, 
                  double idgam[3]);
void opt_lb(double lb[1], Rcpp::NumericVector &expt, Rcpp::NumericVector &a, 
            double idgam[3]);

Rcpp::List em(Rcpp::NumericVector& param2, Rcpp::IntegerVector &xvec, Rcpp::IntegerVector &yvec, 
        Rcpp::IntegerVector &freq, int &n, int &se, int &maxiter, double &tol, int showFlag, int bnb);

#endif
