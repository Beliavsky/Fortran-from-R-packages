// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
# include "header/multiDstack.hpp"
#include <atomic>
using namespace Rcpp;
using namespace RcppParallel;




template<typename valtype, typename indtype>
List mFLSSScomoCpp(
    int len, List vr,
    int _d, int dlst, int dl, int dust, int du,
    int N, NumericVector targetr, NumericVector MEr,
    IntegerVector LBr, IntegerVector UBr, int sizeNeeded, double endTime,
    bool useBisearchInFindBounds)
{
  triM <valtype, indtype> mat;
  std::size_t wordSize = mat.containerWordSize(_d, N, len);
  vec<word> matContainer(wordSize);
  mat.alloc(&matContainer[0], _d, len, N);
  mat.make(&matContainer[0], len, vr);


  valtype ***M = mat.mat;


  vec<indtype> intCtnr;
  vec<valtype> valCtnr;
  vec<mPAT<valtype, indtype> > SK;
  {
    std::size_t stackLen = (unsigned)len + 2;
    unsigned biscaleFactor = (unsigned)std::log2(N + 0.0 - len) + 1;
    intCtnr.resize(stackLen * (stackLen + 1) / 2 * 3 * biscaleFactor);
    valCtnr.resize((3 * (std::size_t)_d + (std::size_t)dl +
      (std::size_t)du) * stackLen * biscaleFactor);
    SK.resize(((unsigned)len + 3) * biscaleFactor);
  }


  mPAT<valtype, indtype> *SKbegin = &SK.front();


  // fill the first SK
  {
    SKbegin->beenUpdated = 0;
    SKbegin->MIN = &valCtnr[0];
    SKbegin->MAX = SKbegin->MIN + dl;
    SKbegin->sumLB = SKbegin->MAX + du;
    SKbegin->sumUB = SKbegin->sumLB + _d;
    SKbegin->sumBresv = SKbegin->sumUB + _d;


    SKbegin->LB = &intCtnr[0];
    SKbegin->UB = SKbegin->LB + len;
    SKbegin->Bresv = SKbegin->UB + len;
    SKbegin->len = len;


    for(indtype i = 0; i < len; ++i)
    {
      SKbegin->LB[i] = LBr[i] - 1;
      SKbegin->UB[i] = UBr[i] - 1;
    }


    // std::copy(targetr.begin(), targetr.end(), SKbegin->target);
    {
      double *vst = &targetr[dlst];
      double *mest = &MEr[dlst];
      for(int i = 0; i < dl; ++i)
      {
        SKbegin->MIN[i] = vst[i] - mest[i];
      }


      vst = &targetr[dust];
      mest = &MEr[dust];
      for(int i = 0; i < du; ++i)
      {
        SKbegin->MAX[i] = vst[i] + mest[i];
      }
    }


    iterSum <valtype, indtype> (SKbegin->sumLB, M[0], SKbegin->LB, len, _d);
    iterSum <valtype, indtype> (SKbegin->sumUB, M[0], SKbegin->UB, len, _d);
  }


  mPAT<valtype,indtype> *SKback = SKbegin + 1; //


  vec<vec<indtype> > rst;
  rst.reserve(sizeNeeded + 6);


  std::atomic<int> totalSize(0);


  // std::size_t SKbackOffset = TTTstack <valtype, indtype> (
  //   len, N,
  //   _d, dlst, dl, dust, du,
  //   M, rst, sizeNeeded,
  //   SKbegin, SKback, useBisearchInFindBounds, totalSize, endTime);
  TTTstack <valtype, indtype> (
      len, N, _d, dlst, dl, dust, du, M, rst, sizeNeeded,
    SKbegin, SKback, useBisearchInFindBounds, totalSize, endTime);


  // Rcpp::List memImage = exportAll<valtype, indtype> (
  //   _d, dlst, dl, dust, du,
  //   len, N, SKbackOffset, matContainer, intCtnr, valCtnr, SK, useBisearchInFindBounds);


  if(rst.size() == 0) return List::create();
  List lis(rst.size());
  for(int i = 0, iend = rst.size(); i < iend; ++i)
  {
    IntegerVector tmp(rst[i].size());
    for(int j = 0, jend = tmp.size(); j < jend; ++j)
    {
      tmp[j] = rst[i][j] + 1;
    }
    lis[i] = tmp;
  }


  return lis;
}




// [[Rcpp::export]]
List z_mFLSSScomo(int len, List vr,
                  int d, int dlst, int dl, int dust, int du,
                  NumericVector targetr,
                  NumericVector MEr, IntegerVector LBr, IntegerVector UBr,
                  int sizeNeededForAll,
                  double duration, bool useFloat, bool useBisearchInFindBounds = 0)
{
  NumericVector tmp = vr[0];
  int vlen = tmp.size();
  double endTime = duration * CLOCKS_PER_SEC + std::clock();
  List result;


  if(std::max(vlen, d) < 127)
  {
    if(useFloat) result = mFLSSScomoCpp <float, signed char> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, useBisearchInFindBounds);
    else result = mFLSSScomoCpp <double, signed char> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, useBisearchInFindBounds);
  }
  else if(std::max(vlen, d) < 32767)
  {
    if(useFloat) result = mFLSSScomoCpp <float, short> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, useBisearchInFindBounds);
    else result = mFLSSScomoCpp <double, short> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, useBisearchInFindBounds);
  }
  else
  {
    if(useFloat) result = mFLSSScomoCpp <float, int> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, useBisearchInFindBounds);
    else result = mFLSSScomoCpp <double, int> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, useBisearchInFindBounds);
  }
  return result;
}
































