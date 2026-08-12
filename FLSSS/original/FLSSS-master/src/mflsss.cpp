// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
# include <chrono>
# include "header/dnyTasking.hpp"
# include "header/mflsssOBJ.hpp"
using namespace Rcpp;
using namespace RcppParallel;
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


// List mFLSSScomoParCpp(
//     int len, NumericMatrix vr,
//     int _d, int dlst, int dl, int dust, int du,
//     int N, NumericVector targetr, NumericVector MEr,
//     IntegerVector LBr, IntegerVector UBr, int sizeNeeded,
//     std::chrono::time_point<std::chrono::steady_clock> endTime,
//     int maxCore, int avgThreadLoad, INT *mask)


// z_mFLSSS(maxCore, len, mV, d, 0, dl, d - du, du, targetMat, ME, LB, UB, solutionNeed, tlimit, useBiSrchInFB)
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
// vr has N rows and _d columns
List mFLSSScpp(
    int len, NumericMatrix vr,
    int _d, int dlst, int dl, int dust, int du, int N,
    NumericMatrix targetMat,
    NumericVector MEr,
    IntegerVector LBr, IntegerVector UBr,
    int sizeNeeded,
    std::chrono::time_point<std::chrono::steady_clock> endTime,
    int maxCore, INT *mask)
{
  // Language("print", vr).eval();
  triM<valtype, indtype> mat;
  std::size_t wordSize = mat.containerWordSize(_d, N, len);
  vec<word> matContainer(wordSize);
  mat.alloc(&matContainer[0], _d, len, N);
  mat.make(&matContainer[0], len, vr);
  valtype ***M = mat.mat;


  // do not give a shit about float vs double any more
  valtype *ME = (valtype*)&*MEr.begin(); // treat is as double
  valtype *targets = (valtype*)&targetMat[0]; // each d elements are a new target
  int Ntasks = targetMat.ncol();


  // indtype commonLB[len], commonUB[len];
  vec<indtype> comB(len * (int)2);
  indtype *commonLB = &*comB.begin();
  indtype *commonUB = commonLB + len;
  for(indtype i = 0; i < len; ++i)
  {
    commonLB[i] = LBr[i] - 1;
    commonUB[i] = UBr[i] - 1;
  }


  // int len, List vr, int _d, int dlst, int dl, int dust, int du, int keyInd, int N,
  shared<valtype, indtype> f(len, N, _d, dlst, dl, dust, du, sizeNeeded, endTime, M, mask);


  vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > mflsssTasks(Ntasks);
  vec<dummyContainers<valtype, indtype, mk, useBiSearch> > dummies(maxCore);


  vec<vec<vec<indtype> > > result(maxCore);
  parMflsssOBJbyCore<valtype, indtype, mk, useBiSearch> (
      mflsssTasks, &f, targets, &dummies[0], ME, commonLB, commonUB, result, Ntasks, maxCore);


  int solutionN = 0;
  {
    for(int i = 0, iend = result.size(); i < iend; ++i)
    {
      solutionN += result[i].size();
    }
  }


  List lis(solutionN);
  for(int k = 0, i = 0, iend = result.size(); i < iend; ++i)
  {
    for(int j = 0, jend = result[i].size(); j < jend; ++j)
    {
      IntegerVector tmp(result[i][j].size());
      for(int u = 0, uend = result[i][j].size(); u < uend; ++u)
      {
        tmp[u] = result[i][j][u] + 1;
      }
      lis[k] = tmp;
      ++k;
    }
  }


  return lis;
}




// [[Rcpp::export]]
List z_mFLSSS(int maxCore, int len, NumericMatrix vr,
              NumericVector maskV,
              int d, int dlst, int dl, int dust, int du,
              NumericMatrix targetMat, NumericVector MEr,
              IntegerVector LBr, IntegerVector UBr,
              int sizeNeed, double duration,
              bool useBiSearch = 0)
{
  int vlen = vr.nrow();
  std::chrono::time_point<std::chrono::steady_clock> endTime =
    std::chrono::steady_clock::now() + std::chrono::seconds(std::size_t(duration));
  List result;
  bool mk = maskV.size() > 0;
  INT *mask = nullptr;
  if (mk) mask = (INT*)&maskV[0];
  
  
  if(std::max(vlen, d) < 127)
  {
         if(mk == 0 and useBiSearch == 0) result = mFLSSScpp<double, signed char, 0, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
    else if(mk == 0 and useBiSearch == 1) result = mFLSSScpp<double, signed char, 0, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
    else if(mk == 1 and useBiSearch == 0) result = mFLSSScpp<INT,    signed char, 1, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
    else if(mk == 1 and useBiSearch == 1) result = mFLSSScpp<INT,    signed char, 1, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
  }
  else if(std::max(vlen, d) < 32767)
  {
         if(mk == 0 and useBiSearch == 0) result = mFLSSScpp<double, short, 0, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
    else if(mk == 0 and useBiSearch == 1) result = mFLSSScpp<double, short, 0, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
    else if(mk == 1 and useBiSearch == 0) result = mFLSSScpp<INT,    short, 1, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
    else if(mk == 1 and useBiSearch == 1) result = mFLSSScpp<INT,    short, 1, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
  }
  else
  {
         if(mk == 0 and useBiSearch == 0) result = mFLSSScpp<double, int, 0, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
    else if(mk == 0 and useBiSearch == 1) result = mFLSSScpp<double, int, 0, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
    else if(mk == 1 and useBiSearch == 0) result = mFLSSScpp<INT,    int, 1, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
    else if(mk == 1 and useBiSearch == 1) result = mFLSSScpp<INT,    int, 1, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, endTime, maxCore, mask);
  }
  return result;
}









