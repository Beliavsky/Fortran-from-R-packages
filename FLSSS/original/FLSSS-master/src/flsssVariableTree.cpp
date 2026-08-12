// # pragma once
# include <Rcpp.h>
// # include <fstream>
// # include "header/PATclass.hpp"
// # include "header/oneDoperation.hpp"
// # include "legacy/singleDIO.hpp"
# include "legacy/macros.hpp"
# include "legacy/singleDstack.hpp"
using namespace Rcpp;



namespace legacy {


template<typename valtype, typename indtype>
List FLSSScpp(
    int len, NumericVector v, double target, double ME,
    IntegerVector LB, IntegerVector UB, int sizeNeed, double durationLimit,
    bool useBisearchInFindBounds)
{
  int N = v.size();
  triMoneD<valtype, indtype> mat;
  std::size_t byteSize = mat.containerByteSize(N, len);
  vec<unsigned char> matContainer(byteSize);
  mat.alloc(&matContainer[0], len, N);
  mat.make(&matContainer[0], len, v);


  valtype **M = mat.mat;


  vec<indtype> intCtnr;
  {
    std::size_t stackLen = len + 2;
    intCtnr.resize(stackLen * (stackLen + 1) / 2 * 3);
  }


  vec<PAT<valtype, indtype> > SK(len + 6);
  PAT<valtype, indtype> *SKbegin = &SK.front();


  // fill the first SK
  {
    SKbegin->LB = &intCtnr[0];
    SKbegin->UB = SKbegin->LB + len;
    SKbegin->UBleftResv = SKbegin->UB + len;
    SKbegin->len = len;


    for(indtype i = 0; i < len; ++i)
    {
      SKbegin->LB[i] = LB[i] - 1;
      SKbegin->UB[i] = UB[i] - 1;
    }


    SKbegin->target = target;


    iterSum <valtype, indtype> (SKbegin->sumLB, M[0], SKbegin->LB, len);
    iterSum <valtype, indtype> (SKbegin->sumUB, M[0], SKbegin->UB, len);
  }


  PAT<valtype, indtype> *SKback = SKbegin + 1; //


  vec<vec<indtype> > rst;
  rst.reserve(sizeNeed + 6);


  // std::size_t SKbackOffset =
    TTTstack<valtype, indtype> (
        len, N, M, ME, rst, sizeNeed, durationLimit,
        SKbegin, SKback, useBisearchInFindBounds);


  // list[[2]][[1]]: len, N, SKbackOffset
  // list[[2]][[2]]: ME
  // list[[2]][[3]]: intCtnr
  // list[[2]][[4]]: triMat
  // list[[1]] stores the results


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


  // Rcpp::List memImage = exportAll<valtype, indtype> (
  //   len, N, SKbackOffset, ME, matContainer, intCtnr, SK, useBisearchInFindBounds);
  // return List::create(Named("solution") = lis, Named("memoryImage") = memImage);


  return lis;
}




/*
template<typename valtype, typename indtype>
List FLSSScatchCpp(List memoryImage, int sizeNeed, double durationLimit)
{
  indtype len, N;
  valtype **M, ME;
  PAT<valtype, indtype> *SKback;
  vec<unsigned char> matContainer;
  vec<indtype> intCtnr;
  vec<PAT<valtype, indtype> > SK;
  bool useBiSearchInFindBounds;


  importAll<valtype, indtype> (memoryImage, len, N, M, SKback, ME, matContainer, intCtnr, SK, useBiSearchInFindBounds);
  PAT<valtype, indtype> *SKbegin = &SK.front();


  vec<vec<indtype> > rst;
  rst.reserve(sizeNeed + 6);


  int SKbackOffset = TTTstack<valtype, indtype> (
    len, N, M, ME, rst, sizeNeed, durationLimit, SKbegin, SKback, useBiSearchInFindBounds);


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


  Rcpp::List memImage = exportAll<valtype, indtype> (
    len, N, SKbackOffset, ME, matContainer, intCtnr, SK, useBiSearchInFindBounds);
  return List::create(Named("solution") = lis, Named("memoryImage") = memImage);
}
*/




}



// [[Rcpp::export]]
List z_FLSSSvariableTree(
    int len, NumericVector v, double target, double ME,
    IntegerVector LB, IntegerVector UB,
    int solutionNeed = 1, double tlimit = 60,
    bool useBiSrchInFB = false, bool useFloat = false)
{
  int N = v.size();
  double durationClock = tlimit * CLOCKS_PER_SEC;
  List result;


  if(N < 127)
  {
    if(useFloat) result = legacy::FLSSScpp <float, signed char> (len, v, target, ME, LB, UB, solutionNeed, durationClock, useBiSrchInFB);
    else result = legacy::FLSSScpp <double, signed char> (len, v, target, ME, LB, UB, solutionNeed, durationClock, useBiSrchInFB);
  }
  else if(N < 32767)
  {
    if(useFloat) result = legacy::FLSSScpp <float, short> (len, v, target, ME, LB, UB, solutionNeed, durationClock, useBiSrchInFB);
    else result = legacy::FLSSScpp <double, short> (len, v, target, ME, LB, UB, solutionNeed, durationClock, useBiSrchInFB);
  }
  else
  {
    if(useFloat) result = legacy::FLSSScpp <float, int> (len, v, target, ME, LB, UB, solutionNeed, durationClock, useBiSrchInFB);
    else result = legacy::FLSSScpp <double, int> (len, v, target, ME, LB, UB, solutionNeed, durationClock, useBiSrchInFB);
  }
  return result;
}




/*
// [[Rcpp::export]]
List z_FLSSSimportMemory(List memImage, int solutionNeed = 1, double tlimit = 60)
{
  IntegerVector types = memImage[memImage.size() - 1]; // types are stored in the last vector
  int &itype = types[0];
  int &ntype = types[1];


  double durationClock = tlimit * CLOCKS_PER_SEC;
  List result;


  if(itype == 1)
  {
    if(ntype == 4) result = FLSSScatchCpp<float, signed char> (memImage, solutionNeed, durationClock);
    else result = FLSSScatchCpp<double, signed char> (memImage, solutionNeed, durationClock);
  }
  else if(itype == 2)
  {
    if(ntype == 4) result = FLSSScatchCpp<float, short> (memImage, solutionNeed, durationClock);
    else result = FLSSScatchCpp<double, short> (memImage, solutionNeed, durationClock);
  }
  else
  {
    if(ntype == 4) result = FLSSScatchCpp<float, int> (memImage, solutionNeed, durationClock);
    else result = FLSSScatchCpp<double, int> (memImage, solutionNeed, durationClock);
  }


  return result;
}
*/


































