// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
using namespace Rcpp;
#include "../arbitraryDimFLSSS/intString2intArray.hpp"
#include "../arbitraryDimFLSSS/mvalFindBound.hpp"
#include <bitset>



// [[Rcpp::export]]
List findBound(StringMatrix &Xval, StringVector &tsum,
               int len, int maxCore = 1,
               IntegerVector lb = IntegerVector(0),
               IntegerVector ub = IntegerVector(0))
{
  vec<uint64_t> Xbit, targetSS, largestSS;
  vec<int> order, colOrder;
  int totalBits = 0;
  bool success = stringMatTo64bitIntMat(
    Xval, tsum, len, maxCore, Xbit, targetSS,
    largestSS, order, colOrder, totalBits);
  vec<uint64_t> &SR = largestSS;// Container recycle.
  if(!success) return List::create();
  int d = targetSS.size(), N = Xbit.size() / d;


  vec<int> LB(len);
  vec<int> UB(len);
  if(lb.size() == 0)
  {
    std::iota(LB.begin(), LB.end(), 0);
    std::iota(UB.begin(), UB.end(), N - len);
  }
  else
  {
    IntegerVector lbz = lb - 1;
    LB.assign(lbz.begin(), lbz.end());
    IntegerVector ubz = ub - 1;
    UB.assign(ubz.begin(), ubz.end());
  }


  TriM Mat(&Xbit[0], d, N, len);


  uint64_t *v = &Xbit[0];
  NumericVector tmptarget(d);
  NumericMatrix tmpv(d, N);
  std::memcpy(&tmpv[0], v, sizeof(double) * N * d);
  std::memcpy(&tmptarget[0], &targetSS[0], sizeof(double) * d);
  tmpv.attr("class") = "integer64";
  tmptarget.attr("class") = "integer64";


  vec<uint64_t> sumLB(d, 0), sumUB(d, 0);
  for(int i = 0; i < len; ++i)
  {
    mvalPlus<int> (&sumLB[0], &sumLB[0], &v[LB[i] * d], d);
    mvalPlus<int> (&sumUB[0], &sumUB[0], &v[UB[i] * d], d);
  }


  // Rcout << "Xbit = \n";
  // Rcout << Xbit[0] << ", " << Xbit[1] << "\n";
  // Rcout << Xbit[2] << ", " << Xbit[3] << "\n";
  // Rcout << Xbit[4] << ", " << Xbit[5] << "\n";
  // Rcout << Xbit[6] << ", " << Xbit[7] << "\n\n";
  //
  //
  //
  // Rcout << "Xbit[0] + Xbit[6] = " << Xbit[0] + Xbit[6] << ", ";
  // Rcout << "Xbit[1] + Xbit[7] = " << Xbit[1] + Xbit[7] << "\n";
  // Rcout << "Xbit[0] + Xbit[2] = " << Xbit[0] + Xbit[2] << ", ";
  // Rcout << "Xbit[1] + Xbit[3] = " << Xbit[1] + Xbit[3] << "\n";
  // Rcout << "targetSS = " << targetSS[0] << ", " << targetSS[1] << "\n";
  // Rcout << "sumLB = " << sumLB[0] << ", " << sumLB[1] << "\n";
  // Rcout << "sumUB = " << sumUB[0] << ", " << sumUB[1] << "\n\n";


  SR.assign(SR.size(), 0);
  int findornot = findBoundCpp<int> (
    len, d, &targetSS[0], &LB[0], &sumLB[0], &UB[0], &sumUB[0], Mat.mat, &SR[0]);


  return List::create(Named("foundOrNot") = findornot,
                      Named("lb") = IntegerVector(LB.begin(), LB.end()) + 1,
                      Named("ub") = IntegerVector(UB.begin(), UB.end()) + 1,
                      Named("v") = transpose(tmpv),
                      Named("targetS") = tmptarget,
                      // Named("triMat") = Rtrimat,
                      Named("colOrder") = colOrder,
                      Named("order") = order);
}








