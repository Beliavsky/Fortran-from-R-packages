

// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#define vec std::vector
#include "../arbitraryDimFLSSS/insertBits.hpp"
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector testInsertBits(List X, IntegerVector NeffectiveBits)
{
  vec<vec<uint64_t> > y(X.size());
  for(int i = 0, iend = y.size(); i < iend; ++i)
  {
    IntegerVector x = X[i];
    y[i].resize(x.size() / 2);
    std::memcpy(&y[i][0], &x[0], sizeof(uint64_t) * y[i].size());
  }
  int totalBits = std::accumulate(
    NeffectiveBits.begin(), NeffectiveBits.end(), 0);
  int finalNint32bit = (totalBits + (sizeof(uint32_t) * 8 - 1)) /
    (sizeof(uint32_t) * 8);
  IntegerVector rst(finalNint32bit, 0);
  insertAllBits((uint64_t*)(&rst[0]), y, &NeffectiveBits[0]);
  return rst;
}


