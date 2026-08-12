// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
#include "../arbitraryDimFLSSS/dnyTasking2.hpp"


// [[Rcpp::export]]
double testTbbPara(NumericVector xv, NumericVector yv,
                   double discount = 10, int Nscaler = 30000,
                   int maxCore = 8)
{
  // CharlieThreadPool tpool(maxCore); // Is a global variable baby.
  // while (true);
  double *x = &xv[0]; // Global shared pointer.
  double *y = &yv[0];
  auto emptyFun = [](std::size_t t)->bool { return false; }; // No early return;


  ParaForEarlyReturn(0, xv.size(), [&](std::size_t i, std::size_t t)->bool
  {
    auto kend = int((1 + x[i] + y[i]) * Nscaler);
    for (int k = 0; k < kend; ++k)
    {
      x[i] = (x[i] + y[i]) * y[i] * discount;
      y[i] = (y[i] + x[i]) * x[i] * discount;
    }
    return false; // No eartly return;
  }, maxCore, 1000, emptyFun, emptyFun);


  std::vector<double> thesum(maxCore, 0);
  double *z = &thesum[0];
  ParaForEarlyReturn(0, xv.size(), [&](std::size_t i, std::size_t t)->bool
  {
    auto kend = int((1 + x[i] + y[i]) * Nscaler);
    for (int k = 0; k < kend; ++k)
    {
      z[t] += (x[i] + y[i]) * y[i] * discount;
      z[t] += (y[i] + x[i]) * x[i] * discount;
    }
    return false; // No eartly return;
  }, maxCore, 1000, emptyFun, emptyFun);


  return std::accumulate(thesum.begin(), thesum.end(), 0.0);


}



