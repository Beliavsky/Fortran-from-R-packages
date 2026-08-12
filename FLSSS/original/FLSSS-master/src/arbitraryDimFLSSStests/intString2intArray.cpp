

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../arbitraryDimFLSSS/intString2intArray.hpp"


// [[Rcpp::export]]
List stringMatTo64bitIntMatTest(
    StringMatrix &Xval, StringVector &tsum, int len, int maxCore)
{
  vec<uint64_t> Xbit, targetSS, largestSS;
  vec<int> order, colOrder;
  int totalBits = 0;
  bool success = stringMatTo64bitIntMat(
    Xval, tsum, len, maxCore, Xbit, targetSS,
    largestSS, order, colOrder, totalBits);
  if(!success) return List::create();
  IntegerMatrix XbitR(largestSS.size() * 2, Xbit.size() / largestSS.size());
  IntegerVector targetSSR(targetSS.size() * 2);
  IntegerVector largestSSR(largestSS.size() * 2);
  std::memcpy(&XbitR[0], &Xbit[0], Xbit.size() * sizeof(uint64_t));
  std::memcpy(&targetSSR[0], &targetSS[0], targetSS.size() * sizeof(uint64_t));
  std::memcpy(&largestSSR[0], &largestSS[0], largestSS.size() * sizeof(uint64_t));
  return List::create(
    Named("NeffectiveBitInRow") = totalBits,
    Named("order") = IntegerVector(order.begin(), order.end()) + 1,
    Named("Xbit") = transpose(XbitR), Named("targetBit") = targetSSR,
    Named("largestSSbit") = largestSSR,
    Named("colOrder") = IntegerVector(colOrder.begin(), colOrder.end()) + 1);
}




