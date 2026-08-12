// [[Rcpp::depends(RcppParallel)]]
# include <fstream>
# include <Rcpp.h>
# include <RcppParallel.h>
// # include "header/multiDstack.hpp"
# include "header/mflsssOBJ.hpp"
using namespace Rcpp;
using namespace RcppParallel;


template<typename valtype, typename indtype, bool mk, bool useBiSearch>
IntegerVector GknapsackCpp(
    int len, NumericMatrix vr, NumericVector profitVec,
    int _d, int dlst, int dl, int dust, int du,
    int N, NumericMatrix targetMat, NumericVector MEr,
    IntegerVector LBr, IntegerVector UBr,
    std::chrono::time_point<std::chrono::steady_clock> endTime, INT *mask,
    int maxCore, int avgThreadLoad, bool verbose, bool approx)
{
  triM<valtype, indtype> mat;
  std::size_t wordSize = mat.containerWordSize(_d, N, len);
  vec<word> matContainer(wordSize);
  mat.alloc(&matContainer[0], _d, len, N);
  mat.make(&matContainer[0], len, vr);


  vec<indtype> optimalSolution(len);
  shared<valtype, indtype> f(
      len, N, _d, dlst, dl, dust, du,
      endTime, mat.mat, mask, // (valtype*)&profitVec[0],
      &profitVec[0],
      &*optimalSolution.begin());


  vec<indtype> acntr(int(2) * len);
  indtype *LBtmp = &*acntr.begin(), *UBtmp = LBtmp + len;
  for(indtype i = 0, iend = len; i < iend; ++i)
  {
    LBtmp[i] = LBr[i] - 1;
    UBtmp[i] = UBr[i] - 1;
  }


  for(int i = 0, iend = targetMat.ncol(); i < iend; ++i)
  {
    valtype *target = (valtype*)&targetMat[0] + i * _d;
    // std::cout << "\n\ntarget = " << target[0] << ", " << target[1] << "\n";

    vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > descendants;
    mitosisForKnapsack<valtype, indtype, mk, useBiSearch> (
        descendants, f, LBtmp, UBtmp, target, (valtype*)&MEr[0],
        maxCore, avgThreadLoad, verbose);


    valtype previousProfit = f.optimalProfit;
    parMflsssOBJforKnapsack<valtype, indtype, mk, useBiSearch> (descendants, maxCore);
    if(f.optimalProfit > previousProfit)
    {
      if(verbose) Rcout << "Updated profit = " << f.optimalProfit << "\n";
      if(approx) break;
    }
  }


  IntegerVector optSolution(f.optimalSolution, f.optimalSolution + len);
  return optSolution + 1;
}




// [[Rcpp::export]]
IntegerVector z_Gknapsack(
  int len, NumericMatrix vr, NumericVector maskV, NumericVector profitVec,
  NumericMatrix targetMat, NumericVector MEr,
  IntegerVector LBr, IntegerVector UBr, double duration,
  bool useBiSearch, int maxCore,
  int avgThreadLoad, bool verbose, bool approx)
{
  int N = profitVec.size();
  int d = vr.ncol();
  std::chrono::time_point<std::chrono::steady_clock> endTime =
    std::chrono::steady_clock::now() + std::chrono::seconds(std::size_t(duration));
  IntegerVector result;
  bool mk = maskV.size() > 0;
  INT* mask = nullptr;
  if (mk) mask = (INT*)&maskV[0];
  if(std::max(N, d) < 127)
  {
         if(mk == 0 and useBiSearch == 0) result = GknapsackCpp<double, signed char, 0, 0> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
    else if(mk == 0 and useBiSearch == 1) result = GknapsackCpp<double, signed char, 0, 1> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
    else if(mk == 1 and useBiSearch == 0) result = GknapsackCpp<INT,    signed char, 1, 0> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
    else if(mk == 1 and useBiSearch == 1) result = GknapsackCpp<INT,    signed char, 1, 1> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
  }
  else if(std::max(N, d) < 32767)
  {
         if(mk == 0 and useBiSearch == 0) result = GknapsackCpp<double, short, 0, 0> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
    else if(mk == 0 and useBiSearch == 1) result = GknapsackCpp<double, short, 0, 1> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
    else if(mk == 1 and useBiSearch == 0) result = GknapsackCpp<INT,    short, 1, 0> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
    else if(mk == 1 and useBiSearch == 1) result = GknapsackCpp<INT,    short, 1, 1> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
  }
  else
  {
         if(mk == 0 and useBiSearch == 0) result = GknapsackCpp<double, int, 0, 0> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
    else if(mk == 0 and useBiSearch == 1) result = GknapsackCpp<double, int, 0, 1> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
    else if(mk == 1 and useBiSearch == 0) result = GknapsackCpp<INT,    int, 1, 0> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
    else if(mk == 1 and useBiSearch == 1) result = GknapsackCpp<INT,    int, 1, 1> (
      len, vr, profitVec, d, 0, 1, 0, d, N, targetMat, MEr, LBr, UBr, endTime, mask,
      maxCore, avgThreadLoad, verbose, approx);
  }
  return result;
}





























