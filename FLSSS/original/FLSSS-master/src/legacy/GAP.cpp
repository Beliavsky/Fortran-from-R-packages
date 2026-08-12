// # pragma once
// [[Rcpp::depends(RcppParallel)]]
# include <Rcpp.h>
# include <RcppParallel.h>
# include <fstream>
// # include <iostream>
// # include <ctime>
// # include <setjmp.h>
# include "header/oneDoperation.hpp"
# include "header/mvalOperation.hpp"
# include "header/macros.hpp"
# include "header/dnyTasking.hpp"
using namespace Rcpp;
// jmp_buf env;


// **v points at this compartment
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
inline unsigned char LBiFind(indtype &lbi, valtype **vst, indtype ubi, valtype *SR,
                             indtype cmpst, indtype dcmp, INT *mask)
{
  if(useBiSearch)
  {
    if(notAllGreaterEqual<valtype, indtype, mk> (
        vst[ubi] + cmpst, SR + cmpst, dcmp, mask)) return 0;
    lbi = mvalLowerBoundBiMan<valtype, indtype, mk> (
      vst + lbi, vst + ubi + 1, SR, cmpst, dcmp, mask) - vst;
  }
  else
  {
    indtype ic = 0;
    lbi = mvalLowerBoundLr<valtype, indtype, mk> (
      vst + lbi, vst + ubi + 1, SR, ic, cmpst, dcmp, mask) - vst;
  }
  if(lbi > ubi) return 0;
  return 1;
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
inline unsigned char UBiFind(indtype &ubi, valtype **vst, indtype lbi, valtype *SR,
                             indtype cmpst, indtype dcmp, INT *mask)
{
  if(useBiSearch)
  {
    if(notAllLessEqual<valtype, indtype, mk> (
        vst[lbi] + cmpst, SR + cmpst, dcmp, mask)) return 0; // Do not delete this
    ubi = mvalUpperBoundBiMan<valtype, indtype, mk> (
      vst + lbi, vst + ubi + 1, SR, cmpst, dcmp, mask) - 1 - vst;
  }
  else
  {
    indtype ic = 0;
    ubi = mvalUpperBoundLr<valtype, indtype, mk> (
      vst + lbi, vst + ubi + 1, SR, ic, cmpst, dcmp, mask) - vst;
  }
  if(ubi < lbi) return 0;
  return 1;
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
inline indtype findBoundCpp(indtype len, indtype d,
                            indtype dlst, indtype dl,
                            indtype dust, indtype du,
                            valtype *Min, valtype *Max,
                            indtype *LB, valtype *sumLB,
                            indtype *UB, valtype *sumUB,
                            valtype **v, INT *mask)
{
  bool boo = false;
  unsigned LBsum = 0, UBsum = 0;
  vec<valtype> acntr((int)2 * d);
  while(true)
  {
    bool boundChanged = false;
    indtype I = 0;
    // valtype cntr[d], *minLessUB = cntr, SR[d];
    valtype *cntr = &*acntr.begin(), *minLessUB = cntr, *SR = cntr + d;
    mvalMinus(minLessUB + dlst, Min, sumUB + dlst, dl);


    LBsum = 0;
    std::fill(sumLB, sumLB + d, 0);
    for(; I < len; ++I)
    {
      indtype LBI = LB[I];
      mvalPlus(SR + dlst, minLessUB + dlst, v[UB[I]] + dlst, dl);
      bool b = LBiFind<valtype, indtype, mk, useBiSearch> (
        LB[I], v, UB[I], SR, dlst, dl, mask);
      if(!b) return 0;
      if(!boundChanged) boundChanged = (LBI != LB[I]);
      mvalPlus(sumLB, sumLB, v[LB[I]], d);
      LBsum += LB[I];
    }


    // Can I stop now
    {
      if(!boo) boo = 1;
      else
      {
        if(!boundChanged)
        {
          if(LBsum == UBsum) // double insurance, also faster
            return 2;
          break;
        }
      }
    }


    I = len - 1;
    // mvalPlusMinus(SR + dust, Max, v[LB[I]] + dust, sumLB + dust, du);
    valtype *maxLessLB = cntr;
    mvalMinus(maxLessLB + dust, Max, sumLB + dust, du);
    boundChanged = 0;


    UBsum = 0;
    std::fill(sumUB, sumUB + d, 0);
    for(; I >= 0; --I)
    {
      indtype UBI = UB[I];
      mvalPlus(SR + dust, maxLessLB + dust, v[LB[I]] + dust, du);
      bool b = UBiFind<valtype, indtype, mk, useBiSearch> (
        UB[I], v, LB[I], SR, dust, du, mask);
      if(!b) return 0;
      if(!boundChanged) boundChanged = (UBI != UB[I]);
      mvalPlus(sumUB, sumUB, v[UB[I]], d);
      UBsum += UB[I];
    }


    // can I jump out now
    {
      if(!boundChanged)
      {
        if(LBsum == UBsum) // double insurance, should be faster too
          return 2;
        break;
      }
    }


  }
  return 1;
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
struct gapPAT
{
  indtype position;
  indtype s, send;
  indtype len;
  indtype *LB, *UB;
  valtype *MIN, *MAX; // Min is of size dl, Max is of size du
  valtype *sumLB, *sumUB;
  double accProfit;


  void copyAnother(gapPAT &X, indtype *XindvStart, valtype *XvalvStart,
                  indtype *indstart, valtype *valstart,
                  indtype d, indtype dl, indtype du)
  {
    position = X.position;
    s = X.s;
    send = X.send;
    len = X.len;
    LB = X.LB - XindvStart + indstart;
    UB = X.UB - XindvStart + indstart;
    MIN = X.MIN - XvalvStart + valstart;
    MAX = X.MAX - XvalvStart + valstart;
    sumLB = X.sumLB - XvalvStart + valstart;
    sumUB = X.sumUB - XvalvStart + valstart;
    std::copy(X.LB, X.LB + len, LB);
    std::copy(X.UB, X.UB + len, UB);
    std::copy(X.MIN, X.MIN + dl + du + d + d, MIN);
    accProfit = X.accProfit;
    // std::memcpy(MIN, X.MIN, sizeof(valtype) * ((std::size_t)d * 2 + dl + du));
  }


  inline void print(indtype d, indtype dl, indtype du)
  {
    Rcpp::Rcout << "position = " << (int)position << ", s = "<< (int)s <<", send = " << (int)send <<
      ", len = " << (int)len << "\n";
    Rcpp::Rcout << "target LB and UB = ";
    for(int i = 0; i < dl; ++i) Rcpp::Rcout <<  MIN[i] << ",";
    Rcpp::Rcout << ",,";
    for(int i = 0; i < du; ++i) Rcpp::Rcout <<  MAX[i] << ", ";
    Rcpp::Rcout << "\n";


    Rcpp::Rcout << "sumLB = ";
    for(int i = 0, iend = d; i < iend; ++i)
    {
      Rcpp::Rcout << sumLB[i] << ", ";
    }
    Rcpp::Rcout << "\n";


    Rcpp::Rcout << "sumUB = ";
    for(int i = 0, iend = d; i < iend; ++i)
    {
      Rcpp::Rcout << sumUB[i] << ", ";
    }
    Rcpp::Rcout << "\n";


    Rcpp::Rcout << "LB = ";
    for(int i = 0, iend = len; i < iend; ++i)
    {
      Rcpp::Rcout << (int)LB[i] << ", ";
    }
    Rcpp::Rcout<<"\n";


    Rcpp::Rcout<<"UB = ";
    for(int i = 0, iend = len; i < iend; ++i)
    {
      Rcpp::Rcout << (int)UB[i] << ", ";
    }
    Rcpp::Rcout << "\n\n";
  }


  inline void print(indtype d, indtype dl, indtype du, std::ofstream &outfile)
  {
    outfile << "position =, " << (int)position << ", s =, "<< (int)s <<", send =, " << (int)send <<
      ", len =, " << (int)len << "\n";
    outfile << "MIN and MAX = ,";
    for(indtype i = 0; i < dl; ++i) outfile <<  MIN[i] << ",";
    outfile << ",,";
    for(indtype i = 0; i < du; ++i) outfile <<  MAX[i] << ", ";
    outfile << "\n";


    outfile << "sumLB =, ";
    for(int i = 0, iend = d; i < iend; ++i)
    {
      outfile << sumLB[i] << ", ";
    }
    outfile << "\n";


    outfile << "sumUB =, ";
    for(int i = 0, iend = d; i < iend; ++i)
    {
      outfile << sumUB[i] << ", ";
    }
    outfile << "\n";


    outfile << "LB =, ";
    for(int i = 0, iend = len; i < iend; ++i)
    {
      outfile << (int)LB[i] << ", ";
    }
    outfile << "\n";


    outfile << "UB =, ";
    for(int i = 0, iend = len; i < iend; ++i)
    {
      outfile << (int)UB[i] << ", ";
    }
    outfile << "\n\n";
  }


  // len is the parent's subset size
  inline void copyParentGene(gapPAT &x, indtype d, indtype dl, indtype du) // x is the parent
  {
    len = x.len;
    MIN = x.sumUB + d;
    MAX = MIN + dl;
    sumLB = MAX + du;
    sumUB = sumLB + d;


    LB = x.UB + len;
    UB = LB + len;


    // std::memcpy(target, x.target, 3 * (std::size_t)d * sizeof(valtype));
    std::memcpy(MIN, x.MIN, sizeof(valtype) * ((std::size_t)dl +
      (std::size_t)du + (std::size_t)d * 2));
    std::memcpy(LB, x.LB, sizeof(indtype) * len); // ! do not think x.LB and x.UB are continous !
    std::memcpy(UB, x.UB, sizeof(indtype) * len);
    accProfit = x.accProfit;
  }


  // equavalent to giveBirth(), and len here is still gene from the parent
  // inline indtype grow(valtype *ME, valtype ***M, indtype d, bool useBiSearch, std::ofstream *outfile = nullptr)
  inline indtype grow(valtype **v, indtype d, indtype dlst, indtype dl,
                      indtype dust, indtype du, INT *mask,
                      double optProfit, double *profitVec,
                      std::ofstream *outfile = nullptr)
  {
    indtype boo = findBoundCpp<valtype, indtype, mk, useBiSearch> (
      len, d, dlst, dl, dust, du, MIN, MAX, LB, sumLB, UB, sumUB, v, mask);


    if(outfile != nullptr)
    {
      *outfile << "Bounds found ___________________________________boo = " << (int)boo << "\n\n";
      print(d, dl, du, *outfile);
    }


    if(boo == 0) return 0;
    if(len == 1) return 3;
    if(boo == 2) return 2;


    if(optProfit > 0)
    {
      double S = accProfit;
      for(indtype i = 0; i < len; ++i)
      {
        S += profitVec[UB[i]];
      }
      if(S <= optProfit) return 0;
    }


    // find the slot that has the least gap
    position = 0;
    indtype Min = *UB - *LB;
    for(indtype i = 1; i < len; ++i)
    {
      indtype tmp = UB[i] - LB[i];
      if(Min > tmp)
      {
        Min = tmp;
        position = i;
      }
    }


    s = UB[position];
    send = LB[position];


    mvalMinus(MIN, MIN, v[s] + dlst, dl);
    mvalMinus(MAX, MAX, v[s] + dust, du);


    mvalMinus(sumLB, sumLB, v[send], d);
    mvalMinus(sumUB, sumUB, v[s], d);


    accProfit += profitVec[s];


    if(position >= len / 2)
    {
      std::copy(LB + position + 1, LB + len, LB + position);
      std::copy(UB + position + 1, UB + len, UB + position);
    }
    else
    {
      std::copy_backward(LB, LB + position, LB + position + 1);
      std::copy_backward(UB, UB + position, UB + position + 1);
      ++LB;
      ++UB;
    }


    --len;
    return 1;
  }


  inline indtype update(valtype **v, indtype d, indtype dlst,
                        indtype dl, indtype dust, indtype du,
                        double *profitVec)
  {
    if(s <= send) return 0;
    --s;
    mvalPlusMinus(MIN, MIN, v[s + 1] + dlst, v[s] + dlst, dl);
    mvalPlusMinus(MAX, MAX, v[s + 1] + dust, v[s] + dust, du);
    accProfit += profitVec[s + 1] - profitVec[s];
    return 1;
  }
};




template<typename valtype, typename indtype>
inline double asolutionProfit(double *profitVec, indtype *asolution, indtype len)
{
  double S = 0;
  for(indtype i = 0; i < len; ++i)
  {
    S += profitVec[asolution[i]];
  }
  return S;
}


// return -1 means time is up
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
signed char TTTstack(
    indtype LEN, indtype N,
    indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
    valtype **V, double *profitVec,
    indtype *optimalRst, // length is LEN
    double *optimalSolProfit,
    gapPAT<valtype, indtype, mk, useBiSearch> *SK,
    gapPAT<valtype, indtype, mk, useBiSearch> *&SKback,
    double endTime, bool verbose, tbb::spin_mutex *mx, INT *mask)
{
  if(SKback <= SK) return 0;


  // std::ofstream outfile("proboutput.csv", std::ios_base::app);
  // outfile << "\n\n\nNew target--------------------------------------------------------\n";
  while(true)
  {
    // std::cout << "rstCurrentSize = " << rstCurrentSize << ", ";
    // (SKback - 1)->print(d, dl, du, outfile);
    // outfile << "parent printed ___________________________________\n\n";


    SKback->copyParentGene(*(SKback - 1), d, dl, du);


    // SKback->print(d, dl, du, outfile);
    // outfile << "parent copied ___________________________________\n\n";


    // indtype boo = SKback->grow(ME, M, d, useBisearchInFindBounds);
    indtype boo = SKback->grow(V, d, dlst, dl, dust, du, mask,
                               *optimalSolProfit, profitVec //, &outfile
                              );


    // outfile << "returned boo =," << (int)boo << "\n-----------";


    // SKback->print(d, dl, du, outfile);
    // outfile << "child grown ___________________________________\n\n";


    // continue to give birth.
    if(boo == 1)
    {
      ++SKback;
      continue;
    }


    // if len in the child becomes 1
    // if(boo == 3 or boo == 2)
    if(boo != 0)
    {
      // indtype common[LEN];
      vec<indtype> acntr(LEN);
      indtype *common = &*acntr.begin();
      for(int i = 1, iend = SKback - &SK[0]; i < iend; ++i)
      {
        common[i - 1] = SK[i].s;
      }
      std::copy(SKback->UB, SKback->UB + SKback->len, common + (SKback - &SK[0] - 1));
      double tmpProfit = asolutionProfit<valtype, indtype> (profitVec, common, LEN);
      mx->lock();
      {
        if(tmpProfit > *optimalSolProfit)
        {
          std::copy(common, common + LEN, optimalRst);
          *optimalSolProfit = tmpProfit;
          // if(verbose) Rcout << "Updated profit = " << tmpProfit << "\n";
        }
      }
      mx->unlock();
    }


    while(true)
    {
      bool updateBool = (SKback - 1)->update(
        V, d, dlst, dl, dust, du, profitVec);


      // (SKback - 1)->print(d, dl, du, outfile);
      // outfile << "parent updated ___________________________________\n\n";


      if(updateBool != 0) break;
      --SKback;
      if(SKback - SK <= 1)
      {
        return 0; // all the combinations have been tried
      }
    }


    if((double)std::clock() > endTime) return -1; // global time up
  }
  return 1;
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
struct parMgap: public RcppParallel::Worker
{
  // bool useBisearchInFindBounds;
  bool verbose;
  indtype len, N;
  indtype d; // d is the dimensionality after padding the key column
  indtype dlst, dl, dust, du;
  valtype **V;
  double *profitVec;
  double endTime;
  INT *mask;
  vec<vec<gapPAT<valtype, indtype, mk, useBiSearch> > > &SKgroup;
  vec<gapPAT<valtype, indtype, mk, useBiSearch>*> &SKbackGroup;
  indtype *optimalSolution;
  double *optimalSolProfit;
  tbb::spin_mutex *mx;
  dynamicTasking *dT;


  void operator() (std::size_t st, std::size_t end)
  {
    // signed char stackRemain;
    for(;;)
    {
      std::size_t objI = 0;
      if(!dT->nextTaskID(objI))
      {
        break;
      }
      // __________________________________________________________________________________________
      // thread function:
      {
        vec<gapPAT<valtype, indtype, mk, useBiSearch> > &SK = SKgroup[objI];
        gapPAT<valtype, indtype, mk, useBiSearch> *SKbegin =
          &SK.front(), *&SKback = SKbackGroup[objI];


        int tmp = TTTstack<valtype, indtype, mk, useBiSearch> (
          len, N, d, dlst, dl, dust, du, V, profitVec,
          optimalSolution, optimalSolProfit, SKbegin, SKback,
          endTime, verbose, mx, mask);
        if(tmp < 0) break;
      }
      // __________________________________________________________________________________________
    }
  }


  parMgap(bool verbose, indtype len, indtype N, indtype d,
          indtype dlst, indtype dl, indtype dust, indtype du,
          double *profitVec, valtype **V, double endTime, INT *mask,
          vec<vec<gapPAT<valtype, indtype, mk, useBiSearch> > > &SKgroup,
          vec<gapPAT<valtype, indtype, mk, useBiSearch>*> &SKbackGroup,
          indtype *optimalSolution, double *optimalSolProfit, int maxCore, int tasks):
    verbose(verbose), len(len), N(N), d(d), dlst(dlst), dl(dl), dust(dust), du(du),
    V(V), profitVec(profitVec), endTime(endTime), mask(mask), SKgroup(SKgroup),
    SKbackGroup(SKbackGroup), optimalSolution(optimalSolution),
    optimalSolProfit(optimalSolProfit) // , outfile(outfile)
  {
    tbb::spin_mutex mxP;
    mx = &mxP;
    dynamicTasking dtask(maxCore, tasks);
    dT = &dtask;
      parallelFor(0, dT->NofCore, *this);
  }
};




/*
template<typename valtype, typename indtype>
// vr has N rows and _d columns
List flsssCompartment(
    int len, // List vr,
    NumericMatrix vr, // each column is an element
    int _d, int dlst, int dl, int dust, int du, int keyInd, int N,
    NumericVector scaleFactorr,
    NumericVector originalTargetr,
    NumericVector keyTargetr,
    NumericVector MEr,
    IntegerVector zeroBasedLB, IntegerVector zeroBasedUB,
    int sizeNeeded, double endTime, int maxCore,
    bool useBiSearchInFindBounds)
{
  vec<valtype*> vptrct(N);
  valtype **V = &vptrct[0];
  vec<valtype> vcontainer;
  if(sizeof(valtype) < 8)
  {
    vcontainer.assign(vr.begin(), vr.end());
    int j = 0;
    for(valtype *i = &vcontainer[0], *iend = &*vcontainer.end(); i < iend; i += _d)
    {
      V[j] = i;
      ++j;
    }
  }
  else
  {
    int j = 0;
    for(valtype *i = (valtype*)&vr[0], *iend = (valtype*)&*vr.end(); i < iend; i += _d)
    {
      V[j] = i;
      ++j;
    }
  }


  valtype *ME = (valtype*)&*MEr.begin();
  vec<valtype> MEcontain;
  if(sizeof(valtype) < 8)
  {
    MEcontain.assign(MEr.begin(), MEr.end());
    ME = &*MEcontain.begin();
  }


  valtype *originalTarget = (valtype*)&*originalTargetr.begin();
  vec<valtype> originalTargetContain;
  if(sizeof(valtype) < 8)
  {
    originalTargetContain.assign(originalTargetr.begin(), originalTargetr.end());
    originalTarget = &*originalTargetContain.begin();
  }


  valtype *keyTarget = (valtype*)&*keyTargetr.begin();
  vec<valtype> keyTargetContain;
  if(sizeof(valtype) < 8)
  {
    keyTargetContain.assign(keyTargetr.begin(), keyTargetr.end());
    keyTarget = &*keyTargetContain.begin();
  }


  valtype *scaleFactor = (valtype*)&*scaleFactorr.begin();
  vec<valtype> scaleFactorContain;
  if(sizeof(valtype) < 8)
  {
    scaleFactorContain.assign(scaleFactorr.begin(), scaleFactorr.end());
    scaleFactor = &*scaleFactorContain.begin();
  }


  // indtype commonLB[len], commonUB[len];
  // for(indtype i = 0; i < len; ++i)
  // {
  //   commonLB[i] = LBr[i] - 1;
  //   commonUB[i] = UBr[i] - 1;
  // }


  vec<vec<vec<indtype> > > result(maxCore);
  vec<vec<indtype> > intCtnrGroup(maxCore);
  vec<vec<valtype> > valCtnrGroup(maxCore);
  vec<vec<gapPAT<valtype, indtype> > > SKgroup(maxCore);
  {
    std::size_t stackLen = len + 2;
    for(int i = 0; i < maxCore; ++i)
    {
      intCtnrGroup[i].resize(stackLen * (stackLen + 1) / 2 * 3);
      valCtnrGroup[i].resize(  stackLen * ((std::size_t)_d * 3 + (std::size_t)dl + (std::size_t)du)  );
      // contain 4 values: MIN, MAX, sumLB, sumUB. MIN and MAX can be less than _d
      SKgroup[i].resize(len + 6);
      result[i].reserve(sizeNeeded);
    }
  }


  vec<unsigned char> keyTargetHasPotential(keyTargetr.size(), 1);
  parMgap<valtype, indtype> (
      useBiSearchInFindBounds, len, N,
      _d, dlst, dl, dust, du,
      sizeNeeded,
      originalTarget, keyTarget, scaleFactor,
      V, ME, &zeroBasedLB[0], &zeroBasedUB[0], result,
      endTime, intCtnrGroup, valCtnrGroup, SKgroup,
      &keyTargetHasPotential[0],
      maxCore, keyTargetr.size());


  // vec<indtype> rst;
  // for(int i = 0, iend = result.size(); i < iend; ++i)
  // {
  //   if(result[i].size() > 0)
  //   {
  //     std::swap(rst, result[i][0]);
  //   }
  // }
  //
  //
  // IntegerVector hope(rst.begin(), rst.end());
  // return hope + 1;


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
List z_flsssCompartment(int maxCore, int len, NumericMatrix V,
                                 int dlst, int dl, int dust, int du, int keyInd,
                                 NumericVector originalTarget, NumericVector keyTarget,
                                 NumericVector scaleFactor, NumericVector MEr,
                                 IntegerVector zeroBasedLB, IntegerVector zeroBasedUB,
                                 int sizeNeed, double duration, bool useFloat, bool useBiSearchInFindBounds)
{
  int N = V.ncol();
  int d = V.nrow();
  double endTime = (double)std::clock() + duration * CLOCKS_PER_SEC;
  List result;
  if(std::max(N, d) < 127)
  {
    if(useFloat) result = flsssCompartment <float, signed char> (len, V, d, dlst, dl, dust, du, keyInd, N, scaleFactor, originalTarget,
       keyTarget, MEr, zeroBasedLB, zeroBasedUB, sizeNeed, endTime, maxCore, useBiSearchInFindBounds);
    else result = flsssCompartment <double, signed char> (len, V, d, dlst, dl, dust, du, keyInd, N, scaleFactor, originalTarget,
       keyTarget, MEr, zeroBasedLB, zeroBasedUB, sizeNeed, endTime, maxCore, useBiSearchInFindBounds);
  }
  else if(std::max(N, d) < 32767)
  {
    if(useFloat) result = flsssCompartment <float, short> (len, V, d, dlst, dl, dust, du, keyInd, N, scaleFactor, originalTarget,
       keyTarget, MEr, zeroBasedLB, zeroBasedUB, sizeNeed, endTime, maxCore, useBiSearchInFindBounds);
    else result = flsssCompartment <double, short> (len, V, d, dlst, dl, dust, du, keyInd, N, scaleFactor, originalTarget,
       keyTarget, MEr, zeroBasedLB, zeroBasedUB, sizeNeed, endTime, maxCore, useBiSearchInFindBounds);
  }
  else
  {
    if(useFloat) result = flsssCompartment <float, int> (len, V, d, dlst, dl, dust, du, keyInd, N, scaleFactor, originalTarget,
       keyTarget, MEr, zeroBasedLB, zeroBasedUB, sizeNeed, endTime, maxCore, useBiSearchInFindBounds);
    else result = flsssCompartment <double, int> (len, V, d, dlst, dl, dust, du, keyInd, N, scaleFactor, originalTarget,
       keyTarget, MEr, zeroBasedLB, zeroBasedUB, sizeNeed, endTime, maxCore, useBiSearchInFindBounds);
  }
  return result;
}
*/




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
void copySKtriplet(vec<gapPAT<valtype, indtype, mk, useBiSearch> > &SK,
                   unsigned depth, vec<indtype> &indvec, vec<valtype> &valvec,
                   vec<gapPAT<valtype, indtype, mk, useBiSearch> > &SKcopy,
                   vec<indtype> &indvecCopy, vec<valtype> &valvecCopy,
                   indtype d, indtype dl, indtype du)
{
  SKcopy.resize(SK.size());
  indvecCopy.resize(indvec.size());
  valvecCopy.resize(valvec.size());
  for(unsigned i = 0; i < depth; ++i)
  {
    SKcopy[i].copyAnother(
        SK[i], &indvec[0], &valvec[0], &indvecCopy[0], &valvecCopy[0], d, dl, du);
  }
}


// spawn from the first level, BFS
template<typename valtype, typename indtype, bool mk, bool useBiSearch>
int spawn(
    vec<gapPAT<valtype, indtype, mk, useBiSearch> > &SK,
    vec<indtype> &indvec,
    vec<valtype> &valvec,
    vec<vec<gapPAT<valtype, indtype, mk, useBiSearch> > > &SKfamily,
    vec<vec<indtype> > &indvecFamily,
    vec<vec<valtype> > &valvecFamily,
    double *profitVec,
    indtype *optimalSolution,
    double &optimalProfit,
    valtype **V,
    indtype LEN, indtype d, indtype dlst,
    indtype dl, indtype dust, indtype du,
    int maxCore, int threadLoad, INT *mask)
{
  int back = 1;
  SKfamily.resize(1); SKfamily[0].swap(SK);
  indvecFamily.resize(1); indvecFamily[0].swap(indvec);
  valvecFamily.resize(1); valvecFamily[0].swap(valvec);


  while(SKfamily.size() > 0
        and SKfamily.size() < unsigned(maxCore * threadLoad)
        and back < LEN)
  {
    vec<vec<vec<gapPAT<valtype, indtype, mk, useBiSearch> > > >
      SKfamilyGroups(SKfamily.size());
    vec<vec<vec<indtype> > > indvecFamilyGroups(SKfamily.size());
    vec<vec<vec<valtype> > > valvecFamilyGroups(SKfamily.size());


    for(int i = 0, iend = SKfamily.size(); i < iend; ++i)
    {
      vec<gapPAT<valtype, indtype, mk, useBiSearch> > &
        SK = SKfamily[i];
      vec<indtype> &indvec = indvecFamily[i];
      vec<valtype> &valvec = valvecFamily[i];


      SK[back].copyParentGene(SK[back - 1], d, dl, du);
      int boo = SK[back].grow(V, d, dlst, dl, dust, du, mask,
                              optimalProfit, profitVec);
      if(boo == 0)
      {
        continue;
      }
      else if(boo == 3 or boo == 2)
      {
        vec<indtype> commonV(LEN);
        indtype *common = &*commonV.begin();
        for(int k = 1; k < back; ++k)
        {
          common[k - 1] = SK[k].s;
        }
        std::copy(SK[back].UB, SK[back].UB + SK[back].len, common + back - 1);
        valtype tmpProfit = asolutionProfit<valtype, indtype> (profitVec, common, LEN);
        if(tmpProfit > optimalProfit)
        {
          std::copy(common, common + LEN, optimalSolution);
          optimalProfit = tmpProfit;
        }
        continue;
      }


      int siblings = SK[back].s - SK[back].send + 1;


      SKfamilyGroups[i].resize(siblings);
      indvecFamilyGroups[i].resize(siblings);
      valvecFamilyGroups[i].resize(siblings);


      SKfamilyGroups[i][0].swap(SK);
      indvecFamilyGroups[i][0].swap(indvec);
      valvecFamilyGroups[i][0].swap(valvec);


      for(int k = 1; k < siblings; ++k)
      {
        copySKtriplet(
          SKfamilyGroups[i][k - 1], back + 1,
          indvecFamilyGroups[i][k - 1],
          valvecFamilyGroups[i][k - 1],
          SKfamilyGroups[i][k],
          indvecFamilyGroups[i][k], valvecFamilyGroups[i][k],
          d, dl, du);
        SKfamilyGroups[i][k][back].update(V, d, dlst, dl, dust, du, profitVec);
      }


      for(int k = 0; k < siblings; ++k)
      {
        SKfamilyGroups[i][k][back].send = SKfamilyGroups[i][k][back].s;
      }
    }


    int spawnTotal = 0;
    for(int i = 0, iend = SKfamilyGroups.size(); i < iend; ++i)
    {
      spawnTotal += SKfamilyGroups[i].size();
    }
    SKfamily.resize(spawnTotal);
    indvecFamily.resize(spawnTotal);
    valvecFamily.resize(spawnTotal);


    int k = 0;
    for(int i = 0, iend = SKfamilyGroups.size(); i < iend; ++i)
    {
      for(int j = 0, jend = SKfamilyGroups[i].size(); j < jend; ++j)
      {
        SKfamilyGroups[i][j].swap(SKfamily[k]);
        indvecFamilyGroups[i][j].swap(indvecFamily[k]);
        valvecFamilyGroups[i][j].swap(valvecFamily[k]);
        ++k;
      }
    }


    ++back;
  }


  return back;
}




template<typename valtype, typename indtype, bool mk, bool useBiSearch>
// vr has N rows and _d columns
IntegerVector GAPcpp(
    int len, // List vr,
    NumericMatrix vr, // each column is an element
    NumericMatrix targetMat, // each column is target
    NumericVector profitVec,
    NumericVector MEr,
    int dlst, int dl, int dust, int du,
    IntegerVector zeroBasedLB, IntegerVector zeroBasedUB,
    double endTime, int maxCore, int threadLoad,
    bool verbose, INT *mask, bool heuristic)
{
  if(maxCore <= 1) threadLoad = 0;
  int N = vr.ncol();
  int _d = vr.nrow();
  vec<valtype*> vptrct(N);
  valtype **V = &vptrct[0];
  {
    int j = 0;
    for(valtype *i = (valtype*)&vr[0], *iend = (valtype*)&*vr.end(); i < iend; i += _d)
    {
      V[j] = i;
      ++j;
    }
  }


  vec<indtype> optimalSolution(len);
  double optimalProfit = 0;
  valtype *ME = (valtype*)&*MEr.begin(); // always double now


  int NofTargets = targetMat.ncol();
  if(verbose) Rcout << "Mining starts. There are " << NofTargets << " tasks:\n";
  for(int i = 0; i < NofTargets; ++i)
  {
    if(verbose) Rcout << i << ", ";
    if((double)std::clock() > endTime) break;


    std::size_t stackLen = len + 3;
    vec<indtype> intCtnr(stackLen * (stackLen + 1) / 2 * 3, 0);
    vec<valtype> valCtnr(stackLen * ((std::size_t)_d * 3 + dl + du), 0.0);
    vec<gapPAT<valtype, indtype, mk, useBiSearch> > SK(len + 6);


    // fill 1st stack slot
    {
      SK[0].MIN = &valCtnr[0];
      SK[0].MAX = SK[0].MIN + dl;
      SK[0].sumLB = SK[0].MAX + du;
      SK[0].sumUB = SK[0].sumLB + _d;
      SK[0].accProfit = 0;


      SK[0].LB = &intCtnr[0];
      SK[0].UB = SK[0].LB + len;
      SK[0].len = len;


      for(indtype i = 0; i < len; ++i)
      {
        SK[0].LB[i] = zeroBasedLB[i];
        SK[0].UB[i] = zeroBasedUB[i];
      }


      {
        valtype *target = (valtype*)&targetMat[0] + i * _d;
        // Derive MIN and MAX
        valtype *vst = target + dlst;
        valtype *me = ME + dlst;
        for(indtype i = 0; i < dl; ++i)
        {
          SK[0].MIN[i] = vst[i] - me[i];
        }
        vst = target + dust;
        me = ME + dust;
        for(indtype i = 0; i < du; ++i)
        {
          SK[0].MAX[i] = vst[i] + me[i];
        }
      }


      iterSum <valtype, indtype> (SK[0].sumLB, V, SK[0].LB, len, _d);
      iterSum <valtype, indtype> (SK[0].sumUB, V, SK[0].UB, len, _d);
    }


    vec<vec<gapPAT<valtype, indtype, mk, useBiSearch> > > SKfamily;
    vec<vec<indtype> > indvecFamily;
    vec<vec<valtype> > valvecFamily;


    // int back = 0;
    double previousProfit = optimalProfit;
    int back = spawn<valtype, indtype, mk, useBiSearch> (
      SK, intCtnr, valCtnr, SKfamily, indvecFamily, valvecFamily,
      &profitVec[0], &optimalSolution[0], optimalProfit,
      V, len, _d, dlst, dl, dust, du, maxCore, threadLoad, mask);


    if(optimalProfit > previousProfit)
    {
      if(verbose) Rcout << "Updated profit = " << optimalProfit << "\n";
      if(heuristic) break;
    }


    if(back >= len) continue;


    vec<gapPAT<valtype, indtype, mk, useBiSearch>*>
      SKfamilyBack(SKfamily.size());
    for(int i = 0, iend = SKfamily.size(); i < iend; ++i)
    {
      SKfamilyBack[i] = &SKfamily[i][0] + back;
    }


    previousProfit = optimalProfit;
    parMgap<valtype, indtype, mk, useBiSearch> (
      verbose, len, N, _d, dlst, dl, dust, du,
      &profitVec[0], V, endTime, mask,
      SKfamily, SKfamilyBack, &optimalSolution[0],
      &optimalProfit, maxCore, SKfamily.size());


    if(optimalProfit > previousProfit)
    {
      if(verbose) Rcout << "Updated profit = " << optimalProfit << "\n";
      if(heuristic) break;
    }
  }


  IntegerVector rst(optimalSolution.begin(), optimalSolution.end());
  return rst;
}




// every column of V is an observation
// [[Rcpp::export]]
IntegerVector z_GAP(int maxCore, int len, NumericMatrix V,
                    NumericVector maskV,
                    int dlst, int dl, int dust, int du,
                    NumericMatrix targetMat,
                    NumericVector profitVec, NumericVector MEr,
                    IntegerVector zeroBasedLB, IntegerVector zeroBasedUB,
                    double duration, bool useBiSearch,
                    int threadLoad = 8, bool verbose = true,
                    bool heuristic = false)
{
  int N = V.ncol();
  int d = V.nrow();
  double endTime = (double)std::clock() + duration * CLOCKS_PER_SEC;
  IntegerVector result;
  bool mk = maskV.size() > 0;
  INT *mask = nullptr;
  if(mk) mask = (INT*)&maskV[0];


  if(std::max(N, d) < 127)
  {
         if(mk == 0 and useBiSearch == 0) result = GAPcpp<double, signed char, 0, 0> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
    else if(mk == 0 and useBiSearch == 1) result = GAPcpp<double, signed char, 0, 1> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
    else if(mk == 1 and useBiSearch == 0) result = GAPcpp<INT,    signed char, 1, 0> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
    else if(mk == 1 and useBiSearch == 1) result = GAPcpp<INT,    signed char, 1, 1> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
  }
  else if(std::max(N, d) < 32767)
  {
         if(mk == 0 and useBiSearch == 0) result = GAPcpp<double, short, 0, 0> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
    else if(mk == 0 and useBiSearch == 1) result = GAPcpp<double, short, 0, 1> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
    else if(mk == 1 and useBiSearch == 0) result = GAPcpp<INT,    short, 1, 0> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
    else if(mk == 1 and useBiSearch == 1) result = GAPcpp<INT,    short, 1, 1> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
  }
  else
  {
         if(mk == 0 and useBiSearch == 0) result = GAPcpp<double, int, 0, 0> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
    else if(mk == 0 and useBiSearch == 1) result = GAPcpp<double, int, 0, 1> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
    else if(mk == 1 and useBiSearch == 0) result = GAPcpp<INT,    int, 1, 0> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
    else if(mk == 1 and useBiSearch == 1) result = GAPcpp<INT,    int, 1, 1> (
      len, V, targetMat, profitVec, MEr, dlst, dl, dust, du,
      zeroBasedLB, zeroBasedUB, endTime, maxCore, threadLoad, verbose, mask, heuristic);
  }
  return result;
}













































