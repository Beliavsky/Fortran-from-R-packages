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


template<typename T> // x gets destroyed. T must have .swap() method.
void pushbackVdestruct(vec<T> &v, T &x)
{
  if(v.size() == v.capacity()) v.resize(v.size() * 2);
  v.resize(v.size() + 1);
  v.back().swap(x);
}
template<typename T> // x gets destroyed. T must have .swap() method.
void concatenateVdestruct(vec<T> &v, vec<T> &x)
{
  int vsize = v.size();
  std::size_t desiredSize = v.size() + x.size();
  if(desiredSize > v.capacity()) v.resize(desiredSize * 2);
  v.resize(desiredSize);
  for(int i = 0, iend = x.size(); i < iend; ++i)
    v[vsize + i].swap(x[i]);
}
template<typename T> // x gets destroyed. v is the unpacked version of x.
void unpack(vec<T> &v, vec<vec<T> > &x)
{
  int size = 0;
  for(int i = 0, iend = x.size(); i < iend; ++i) size += x[i].size();
  v.resize(size);
  for(int k = 0, i = 0, iend = x.size(); i < iend; ++i)
  {
    for(int j = 0, jend = x[i].size(); j < jend; ++j)
    {
      v[k].swap(x[i][j]);
      ++k;
    }
  }
}




// z_mFLSSS(maxCore, len, mV, d, 0, dl, d - du, du, targetMat, ME, LB, UB, solutionNeed, tlimit, useBiSrchInFB)
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
// vr has N rows and _d columns
List mFLSSScppImage(
    int len, NumericMatrix vr,
    int _d, int dlst, int dl, int dust, int du, int N,
    NumericMatrix targetMat,
    NumericVector MEr,
    IntegerVector LBr, IntegerVector UBr,
    int sizeNeeded,
    INT *mask, int Ninstance)
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


  shared<valtype, indtype> f(len, N, _d, dlst, dl, dust, du, sizeNeeded,
                             std::chrono::steady_clock::now(), M, mask);


  vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > mflsssOBJvec;
  vec<vec<indtype> > result;
  if(Ntasks > Ninstance)
  {
    mflsssOBJvec.resize(Ntasks);
    for(int i = 0, iend = Ntasks; i < iend; ++i)
    {
      mflsssOBJvec[i].initialize(
          &f, targets + i * f.d, ME, commonLB, commonUB);
    }
  }
  else
  {
    mflsssOBJvec.reserve(Ninstance * 2);
    result.reserve(Ninstance * 2);
    // Rcout << "Ninstance = " << Ninstance << "\n";
    // Rcout << "Ntasks = " << Ntasks << "\n";
    // Rcout << "Ninstance / Ntasks + 1 = " << Ninstance / Ntasks + 1 << "\n";
    for(int i = 0, iend = Ntasks; i < iend; ++i)
    {
      vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > descendants;
      vec<vec<indtype> > rst;
      mitosis<valtype, indtype, mk, useBiSearch> (
          descendants, f, rst, commonLB, commonUB, targets + i * f.d,
          ME, Ninstance / Ntasks + 1, 1);
      concatenateVdestruct(mflsssOBJvec, descendants);
      // Rcout << "loop 1, rst.size() = " << rst.size() << "\n";
      concatenateVdestruct(result, rst);
      // Rcout << "loop 1, result.size() = " << result.size() << "\n";
    }
    // Rcout << "mflsssOBJvec.size() = " << mflsssOBJvec.size() << "\n";
    // return List::create();


    while(result.size() < (unsigned)(sizeNeeded) and
          mflsssOBJvec.size() < unsigned(Ninstance) and
          mflsssOBJvec.size() != 0)
    {
      // Rcout << "1.2\n";
      vec<vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > > ds(mflsssOBJvec.size());
      vec<vec<vec<indtype> > > rt; rt.reserve(7);
      // Rcout << "mflsssOBJvec.size() = " << mflsssOBJvec.size() << "\n";
      for(int i = 0, iend = mflsssOBJvec.size(); i < iend; ++i)
      {
        vec<vec<indtype> > rst;
        // Rcout << "i = " << i << "\n";
        mitosis(ds[i], mflsssOBJvec[i], f, rst); // Grow pairs.
        // Rcout << "loop 2, rst.size() = " << rst.size() << "\n";
        pushbackVdestruct(rt, rst);
        // Rcout << "loop 2, rt.size() = " << rt.size() << "\n";
      }
      // Rcout << "1.3\n";


      vec<mflsssOBJ<valtype, indtype, mk, useBiSearch> > mflsssOBJvecTmp;
      unpack(mflsssOBJvecTmp, ds);
      mflsssOBJvecTmp.swap(mflsssOBJvec);


      vec<vec<indtype> > resultTmp;
      unpack(resultTmp, rt);
      concatenateVdestruct(result, resultTmp);
    }
  }


  if constexpr (false)
  {


    Rcout << "f.totalSize < f.sizeNeed = " << int(f.totalSize < f.sizeNeed) << "\n";
    Rcout << "mflsssOBJvec.size() = " << int(mflsssOBJvec.size()) << "\n";
    if(f.totalSize < f.sizeNeed and mflsssOBJvec.size() > 0)
    {
      f.endTime = std::chrono::steady_clock::now() +
        std::chrono::seconds(std::size_t(1e6));
      parMflsssOBJ<valtype, indtype, mk, useBiSearch> (mflsssOBJvec, 15);
    }


    // vec<vec<indtype> > rst; rst.reserve(7);
    for(int i = 0, iend = mflsssOBJvec.size(); i < iend; ++i)
    {
      vec<vec<indtype> > &resulttmp = mflsssOBJvec[i].result;
      // Rcout << "resulttmp.size() = " << int(resulttmp.size()) << "\n";
      for(int j = 0, jend = resulttmp.size(); j < jend; ++j)
      {
        result.push_back(resulttmp[j]);
      }
    }


    List lis(result.size());
    for(int i = 0, iend = result.size(); i < iend; ++i)
    {
      IntegerVector tmp(result[i].begin(), result[i].end());
      lis[i] = tmp + 1;
    }
    return lis;


  }


  if(false)
  {
    for(int i = 0, iend = mflsssOBJvec.size(); i < iend; ++i)
      mflsssOBJvec[i].print();
  }




  // Rcout << "result.size() = " << result.size() << "\n";
  if(mflsssOBJvec.size() == 0)
  {
    List existingRst(result.size());
    for(int i = 0, iend = result.size(); i < iend; ++i)
    {
      IntegerVector tmp(result[i].begin(), result[i].end());
      existingRst[i] = tmp + 1;
    }
    return List::create(
      Named("mflsssObjects") =
        List::create(List::create(Named("mflsssInner") = R_NilValue,
                     Named("sharedSave") = R_NilValue,
                     Named("vr") = vr,
                     Named("indtype") = int(sizeof(indtype)),
                     Named("mk") = mk,
                     Named("useBiSearch") = useBiSearch)),
      Named("solutionsFound") = existingRst);
  }


  f.totalSize = result.size();
  RawVector sharedSave = f.save();
  // vr
  List mflsssOBJlist(mflsssOBJvec.size());
  List existingRst(result.size());
  for(int i = 0, iend = existingRst.size(); i < iend; ++i)
  {
    IntegerVector tmp(result[i].begin(), result[i].end());
    existingRst[i] = tmp + 1;
  }


  for(int i = 0, iend = mflsssOBJlist.size(); i < iend; ++i)
  {
    List lis = mflsssOBJvec[i].save();
    mflsssOBJlist[i] = List::create(Named("mflsssInner") = lis,
                                    Named("sharedSave") = sharedSave,
                                    Named("vr") = vr,
                                    Named("indtype") = int(sizeof(indtype)),
                                    Named("mk") = mk,
                                    Named("useBiSearch") = useBiSearch);
  }
  return List::create(Named("mflsssObjects") = mflsssOBJlist,
                      Named("solutionsFound") = existingRst);
}




// [[Rcpp::export]]
List z_mFLSSSimage(int len, NumericMatrix vr,
                   NumericVector maskV,
                   int d, int dlst, int dl, int dust, int du,
                   NumericMatrix targetMat, NumericVector MEr,
                   IntegerVector LBr, IntegerVector UBr,
                   int sizeNeed, bool useBiSearch = 0,
                   int Ninstance = 100000)
{
  int vlen = vr.nrow();
  List result;
  INT *mask = nullptr;
  bool mk = maskV.size() > 0;
  if (mk) mask = (INT*)&maskV[0];
  
  
  if(std::max(vlen, d) < 127)
  {
    if(mk == 0 and useBiSearch == 0) result = mFLSSScppImage<double, signed char, 0, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
    else if(mk == 0 and useBiSearch == 1) result = mFLSSScppImage<double, signed char, 0, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
    else if(mk == 1 and useBiSearch == 0) result = mFLSSScppImage<INT,    signed char, 1, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
    else if(mk == 1 and useBiSearch == 1) result = mFLSSScppImage<INT,    signed char, 1, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
  }
  else if(std::max(vlen, d) < 32767)
  {
    if(mk == 0 and useBiSearch == 0) result = mFLSSScppImage<double, short, 0, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
    else if(mk == 0 and useBiSearch == 1) result = mFLSSScppImage<double, short, 0, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
    else if(mk == 1 and useBiSearch == 0) result = mFLSSScppImage<INT,    short, 1, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
    else if(mk == 1 and useBiSearch == 1) result = mFLSSScppImage<INT,    short, 1, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
  }
  else
  {
    if(mk == 0 and useBiSearch == 0) result = mFLSSScppImage<double, int, 0, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
    else if(mk == 0 and useBiSearch == 1) result = mFLSSScppImage<double, int, 0, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
    else if(mk == 1 and useBiSearch == 0) result = mFLSSScppImage<INT,    int, 1, 0> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
    else if(mk == 1 and useBiSearch == 1) result = mFLSSScppImage<INT,    int, 1, 1> (
      len, vr, d, dlst, dl, dust, du, vlen, targetMat, MEr, LBr, UBr, sizeNeed, mask, Ninstance);
  }
  return result;
}












// valtype is always double.
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
List mFLSSSimport(List mflsssObj, int sizeNeed, double tlimit)
{
  List mflsssInner = mflsssObj["mflsssInner"];
  RawVector sharedSave = mflsssObj["sharedSave"];
  shared<valtype, indtype> f;
  f.read(sharedSave);
  triM<valtype, indtype> mat;
  std::size_t wordSize = mat.containerWordSize(f.d, f.N, f.subsetSize);
  vec<word> matContainer(wordSize);
  mat.alloc(&matContainer[0], f.d, f.subsetSize, f.N);
  NumericMatrix vr = mflsssObj["vr"];
  mat.make(&matContainer[0], f.subsetSize, vr);
  valtype ***M = mat.mat;
  NumericVector maskV = mflsssObj["maskV"];
  INT* maskVptr = nullptr;
  if (maskV.size() != 0) maskVptr = (INT*)(&maskV[0]);
  f.read(M, maskVptr, nullptr, nullptr);
  f.endTime = std::chrono::steady_clock::now() + std::chrono::seconds(std::size_t(tlimit));
  f.sizeNeed = sizeNeed;
  f.totalSize = 0;


  mflsssOBJ<valtype, indtype, mk, useBiSearch> mob;
  // Rcout << "1.0\n";
  mob.read(mflsssInner, &f);
  // Rcout << "2.0\n";


  // return List::create();
  if(false)
  {
    mob.print();
  }


  mob.TTTstackRun();


  List lis(mob.result.size());
  for(int i = 0, iend = mob.result.size(); i < iend; ++i)
  {
    IntegerVector v(mob.result[i].begin(), mob.result[i].end());
    lis[i] = v + 1;
  }
  return lis;
}


template<typename valtype, typename indtype, bool mk>
List mFLSSSimport(List mflsssObj, int sizeNeed, double tlimit)
{
  int useBiSearch = mflsssObj["useBiSearch"];
  if(useBiSearch) return mFLSSSimport<valtype, indtype, mk, true> (
      mflsssObj, sizeNeed, tlimit);
  return mFLSSSimport<valtype, indtype, mk, false> (
      mflsssObj, sizeNeed, tlimit);
}


template<typename valtype, typename indtype>
List mFLSSSimport(List mflsssObj, int sizeNeed, double tlimit)
{
  int mk = mflsssObj["mk"];
  if(mk) return mFLSSSimport<valtype, indtype, true> (
      mflsssObj, sizeNeed, tlimit);
  return mFLSSSimport<valtype, indtype, false> (
      mflsssObj, sizeNeed, tlimit);
}


List mFLSSSimport(List mflsssObj, int sizeNeed, double tlimit)
{
  List rst;
  int a = mflsssObj["indtype"];
  if(a == 1) rst = mFLSSSimport<double, signed char> (mflsssObj, sizeNeed, tlimit);
  else if(a == 2) rst = mFLSSSimport<double, short> (mflsssObj, sizeNeed, tlimit);
  else if(a == 4) rst = mFLSSSimport<double, int> (mflsssObj, sizeNeed, tlimit);
  return rst;
}


// [[Rcpp::export]]
List z_mFLSSSimport(List mflsssObj, int sizeNeed, double tlimit)
{
  if(mflsssObj.size() == 0) return List::create();
  List mflsssInner = mflsssObj["mflsssInner"];
  if(mflsssInner.size() <= 1) return List::create();
  return mFLSSSimport(mflsssObj, sizeNeed, tlimit);
}








