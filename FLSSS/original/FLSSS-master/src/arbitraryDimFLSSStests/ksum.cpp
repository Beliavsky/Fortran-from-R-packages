// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#define vec std::vector
// #include "../arbitraryDimFLSSS/intString2intArray.hpp"
#include "../arbitraryDimFLSSS/ksumHash.hpp"
using namespace Rcpp;


#ifndef vec
#define vec std::vector
#endif


// Every column of v is an observation.
// [[Rcpp::export]]
RawVector ksumHash(NumericMatrix v, IntegerVector lb, IntegerVector ub,
                   int maxCore, int upscale, bool verbose = true)
{
  vec<uint64_t> x(v.begin(), v.end());
  int k = lb.size(), N = v.ncol(), d = v.nrow();
  Ksum K;
  K.reset(&x[0], N, d, k, &lb[0], &ub[0], upscale, maxCore);
  auto rst = K(verbose);
  return rst;
}


// [[Rcpp::export]]
void checkAddress(RawVector x)
{
  Rcout << uint64_t(&x[10]) << "\n";
}











