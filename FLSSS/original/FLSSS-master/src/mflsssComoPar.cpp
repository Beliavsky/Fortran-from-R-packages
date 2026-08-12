// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
# include <fstream>
# include <Rcpp.h>
# include <RcppParallel.h>
# include "header/mflsssOBJ.hpp"
using namespace Rcpp;
using namespace RcppParallel;
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
List mFLSSScomoParCpp(
    int len, NumericMatrix vr,
    int _d, int dlst, int dl, int dust, int du,
    int N, NumericVector targetr, NumericVector MEr,
    IntegerVector LBr, IntegerVector UBr, int sizeNeeded,
    std::chrono::time_point<std::chrono::steady_clock> endTime,
    int maxCore, int avgThreadLoad, INT *mask)
{
  triM <valtype, indtype> mat;
  std::size_t wordSize = mat.containerWordSize(_d, N, len);
  vec<word> matContainer(wordSize);
  mat.alloc(&matContainer[0], _d, len, N);
  mat.make(&matContainer[0], len, vr);


  shared<valtype, indtype> f(
      len, N, _d, dlst, dl, dust, du, sizeNeeded, endTime, mat.mat, mask);


  // indtype LBtmp[len], UBtmp[len];
  vec<indtype> lbub(len * (int)2);
  indtype *LBtmp = &*lbub.begin(), *UBtmp = LBtmp + len;
  for(indtype i = 0, iend = len; i < iend; ++i)
  {
    LBtmp[i] = LBr[i] - 1;
    UBtmp[i] = UBr[i] - 1;
  }
  // vec<valtype> targetTmp(targetr.begin(), targetr.end());
  // vec<valtype> meTmp(MEr.begin(), MEr.end());


  vec<vec<indtype> > rst;
  rst.reserve(7);


  // Rcout << "1.1\n";
  vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > descendants;
  mitosis<valtype, indtype, mk, useBiSearch> (
      descendants, f, rst, LBtmp, UBtmp, (valtype*)&*targetr.begin(),
      (valtype*)&*MEr.begin(), maxCore, avgThreadLoad);


  // for(int i = 0, iend = rst.size(); i < iend; ++i)
  // {
  //   for(int j = 0, jend = rst[i].size(); j < jend; ++j)
  //   {
  //     Rcout << (int)(rst[i][j]) << ", ";
  //   }
  //   Rcout << "\n";
  // }
  // if(true)
  // {
  //   for(int i = 0, iend = descendants.size(); i < iend; ++i)
  //     descendants[i].print();
  // }


  // Rcout << "descendants.size() = " << descendants.size() << "\n";
  if(f.totalSize < f.sizeNeed and descendants.size() > 0)
  {
    // Rcout << "parallel mflsss object computing\n";
    // Rcout << "N of tasks = " << descendants.size() << "\n";
    parMflsssOBJ<valtype, indtype, mk, useBiSearch> (descendants, maxCore);
  }


  // Rcout << "-1\n";
  for(int i = 0, iend = descendants.size(); i < iend; ++i)
  {
    vec<vec<indtype> > &result = descendants[i].result;
    for(int j = 0, jend = result.size(); j < jend; ++j)
    {
      // for(int k = 0, kend = result[j].size(); k < kend; ++k)
      // {
      //   Rcout << (int)result[j][k] << ", ";
      // }
      // Rcout << "\n";
      rst.push_back(result[j]);
    }
    // Rcout << "-1\n";
  }


  List lis(rst.size());
  for(int i = 0, iend = rst.size(); i < iend; ++i)
  {
    IntegerVector tmp(rst[i].begin(), rst[i].end());
    lis[i] = tmp + 1;
  }
  return lis;
}




// [[Rcpp::export]]
List z_mFLSSScomoPar(int maxCore, int len, NumericMatrix vr, NumericVector maskV,
                  int d, int dlst, int dl, int dust, int du,
                  NumericVector targetr,
                  NumericVector MEr, IntegerVector LBr, IntegerVector UBr,
                  int sizeNeededForAll,
                  double duration, bool useBiSearch = 0, int avgThreadLoad = 8)
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
         if(mk == 0 and useBiSearch == 0) result = mFLSSScomoParCpp <double, signed char, 0, 0> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
    else if(mk == 0 and useBiSearch == 1) result = mFLSSScomoParCpp <double, signed char, 0, 1> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
    else if(mk == 1 and useBiSearch == 0) result = mFLSSScomoParCpp <INT   , signed char, 1, 0> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
    else if(mk == 1 and useBiSearch == 1) result = mFLSSScomoParCpp <INT   , signed char, 1, 1> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
  }
  else if(std::max(vlen, d) < 32767)
  {
         if(mk == 0 and useBiSearch == 0) result = mFLSSScomoParCpp <double, short, 0, 0> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
    else if(mk == 0 and useBiSearch == 1) result = mFLSSScomoParCpp <double, short, 0, 1> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
    else if(mk == 1 and useBiSearch == 0) result = mFLSSScomoParCpp <INT   , short, 1, 0> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
    else if(mk == 1 and useBiSearch == 1) result = mFLSSScomoParCpp <INT   , short, 1, 1> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
  }
  else
  {
         if(mk == 0 and useBiSearch == 0) result = mFLSSScomoParCpp <double, int, 0, 0> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
    else if(mk == 0 and useBiSearch == 1) result = mFLSSScomoParCpp <double, int, 0, 1> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
    else if(mk == 1 and useBiSearch == 0) result = mFLSSScomoParCpp <INT   , int, 1, 0> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
    else if(mk == 1 and useBiSearch == 1) result = mFLSSScomoParCpp <INT   , int, 1, 1> (len, vr, d, dlst, dl, dust, du, vlen, targetr, MEr, LBr, UBr, sizeNeededForAll, endTime, maxCore, avgThreadLoad, mask);
  }
  return result;
}





















































