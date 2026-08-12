#pragma once
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <chrono>
// #include "dnyTasking2.hpp"
#include "charlieThreadPool.hpp"
#include "mflsssOBJ.hpp"
using namespace Rcpp;
#include "intString2intArray.hpp"


template<typename T> // x gets destroyed. T must have .swap() method.
inline void concatenateVdestruct(vec<T> &v, vec<T> &x)
{
  int vsize = v.size();
  std::size_t desiredSize = v.size() + x.size();
  if (desiredSize > v.capacity()) v.resize(desiredSize * 2);
  v.resize(desiredSize);
  for (int i = 0, iend = x.size(); i < iend; ++i)
    v[vsize + i].swap(x[i]);
}


template<typename ing>
inline bool makeMflsssObjs(int len, StringMatrix Xval, StringVector tsum,
                           int sizeNeeded, int Ninstance, CharlieThreadPool &tp,
                           vec<mflsssOBJ<ing> > &mflsssOBJvec,
                           vec<vec<ing> > &result, Shared<ing> &f,
                           double tlimit, vec<int> &order, TriM &Mat,
                           KsumLookUpTable<ing> &ksumtable,
                           bool ksumIsComputedInternallyAndWillNotBeExported)
  // ksumtable is here merely for resetting f.
{


  vec<uint64_t> Xbit, targetSS, largestSS;
  vec<int> colOrder;
  int totalBits = 0;


  bool success = stringMatTo64bitIntMat(
    Xval, tsum, len, tp, Xbit, targetSS, largestSS,
    order, colOrder, totalBits, ksumIsComputedInternallyAndWillNotBeExported);


  if (!success) return false;
  ing d = targetSS.size(), N = Xbit.size() / d;
  Mat.reset(&Xbit[0], d, N, len);


  uint64_t ***M = Mat.mat;
  uint64_t *targets = &targetSS[0];


  vec<ing> comB(len * (int)2);
  ing *commonLB = &*comB.begin();
  ing *commonUB = commonLB + len;
  for(ing i = 0; i < len; ++i)
  {
    commonLB[i] = i;
    commonUB[i] = N - len + i;
  }


  TimePoint endTime = std::chrono::steady_clock::now() +
    std::chrono::seconds(std::size_t(tlimit));
  f.reset(len, N, d, sizeNeeded, endTime, M, &ksumtable);


  Mitosis<ing> mitosis;
  {
    vec<mflsssOBJ<ing> > descendants;
    vec<vec<ing> > rst;
    mitosis(descendants, f, rst, commonLB, commonUB, targets, Ninstance);
    concatenateVdestruct(mflsssOBJvec, descendants);
    concatenateVdestruct(result, rst);
  }


  return true;
}




template<typename ing>
inline void extractBoundsInfoFromMflsssObj(
    mflsssOBJ<ing> *mflsssOBJ, int NflsssObj,
    vec<ing*> &existingLBs,
    vec<ing*> &existingUBs,
    vec<ing>  &existingBsizes)
{
  existingLBs.resize(NflsssObj);
  existingUBs.resize(NflsssObj);
  existingBsizes.resize(NflsssObj);
  for (int i = 0; i < NflsssObj; ++i)
  {
    // Rcpp::Rcout << mflsssOBJ[i].SKback << ", ";
    auto p = mflsssOBJ[i].SKback->parent;
    existingBsizes[i] = p->len;
    existingLBs[i] = p->LB;
    existingUBs[i] = p->UB;
  }
}




// givenKsumTable could be changed.
// Even though it is an Rcpp::List, here we still need to use reference, because
// the object itself is subject to change.
template<typename ing>
inline void processGivenKsumtableORrecompute(
    mflsssOBJ<ing> *X, int NflsssObj, Rcpp::List &givenKsumTable, bool verbose,
    int ksum, int upscale, CharlieThreadPool &tp, KsumLookUpTable<ing> &ksumtable,
    Shared<ing> &f) // f.M should have been initialized.
{
  if (givenKsumTable.size() != 0)
  {
    if (verbose) Rcpp::Rcout << "A k-sum accelerator is given. Read in..\n\n";
    ksumtable.read(givenKsumTable, f.subsetSize);
  }
  else if (ksum >= 3)
  {
    if (verbose) Rcpp::Rcout << "No k-sum accelerator is given. \
User asks to compute one..\n\n";
    vec<ing*> existingLBs, existingUBs;
    vec<ing>  existingBsizes;
    // Rcpp::Rcout << "NflsssObj = " << NflsssObj << "\n";
    extractBoundsInfoFromMflsssObj<ing> (
        X, NflsssObj, existingLBs, existingUBs, existingBsizes);
    givenKsumTable = ksumtable.make(
      existingLBs, existingUBs, existingBsizes, f.subsetSize, f.N, f.d,
      f.M[0][0], ksum, upscale, tp, verbose);
    // if (givenKsumTable.size() == 0) f.ksumtable = nullptr;
    // else f.ksumtable = &ksumtable;
  }
  else
  {
    if (verbose) Rcpp::Rcout << "No k-sum accelerator is given. " <<
      int(ksum) << "-sum accelerator is ignored.\n\n";
    // f.ksumtable = nullptr;
  }
  f.ksumtable = &ksumtable;
}




template<typename ing>
inline List saveMflsssObjs(vec<mflsssOBJ<ing> > &X, Shared<ing> &f,
                           TriM &Mat, vec<int> &order,
                           Rcpp::List ksumtable)
{
  auto fsave = f.save();
  auto ftrimat = Mat.save();
  IntegerVector odr(order.begin(), order.end());
  List rst(X.size());
  for (int i = 0, iend = rst.size(); i < iend; ++i)
  {
    rst[i] = List::create(
      Named("indtype") = int(sizeof(ing)),
      Named("shared") = fsave,
      Named("trimat") = ftrimat,
      Named("order") = odr,  Named("obj") = X[i].save(),
      Named("ksumtable") = ksumtable);
  }
  return rst;
}




template<typename ing>
inline List decomposeArbFLSSSexport(
    int len, StringMatrix V, StringVector target,
    int approxNinstance, CharlieThreadPool &tp,
    SEXP ksumtableList, int ksum, int upscale, bool verbose)
{
  vec<mflsssOBJ<ing> > mflsssOBJvec;
  vec<vec<ing> > result;
  Shared<ing> f;
  vec<int> order;
  TriM Mat;
  KsumLookUpTable<ing> ksumtable(len);


  bool ksumWillBeComputedInternally = TYPEOF(ksumtableList) == NILSXP;


  bool success = makeMflsssObjs<ing> (
    len, V, target, 2147483647, approxNinstance, tp,
    mflsssOBJvec, result, f, 9000000000.0, order, Mat, ksumtable,
    ksumWillBeComputedInternally);
  if (!success) return List();


  List ksumtableConcrete;
  if (!ksumWillBeComputedInternally) ksumtableConcrete = ksumtableList;
  processGivenKsumtableORrecompute<ing> (
    &mflsssOBJvec[0], mflsssOBJvec.size(), ksumtableConcrete, verbose,
    ksum, upscale, tp, ksumtable, f);


  List rstFound(result.size());
  for (int i = 0, iend = result.size(); i < iend; ++i)
  {
    IntegerVector tmp(result[i].size());
    for (int k = 0, kend = tmp.size(); k < kend; ++k)
      tmp[k] = order[result[i][k]] + 1;
    rstFound[i] = tmp;
  }


  return List::create(
    Named("arbFLSSSobjects") = saveMflsssObjs(
      mflsssOBJvec, f, Mat, order, ksumtableConcrete),
    Named("solutionsFound") = rstFound);
}




// Ndecomposed is the number of mflsss objects.
template<typename ing>
inline List mflsssArb(int len, StringMatrix Xval, StringVector tsum,
                      int sizeNeeded, CharlieThreadPool &tp, int Ndecomposed,
                      double tlimit, List givenKsumTable, // Can outsource the table.
                      int ksum, int upscale, bool verbose)
{


  vec<mflsssOBJ<ing> > mflsssOBJvec;
  TriM Mat;
  vec<vec<ing> > mitosisRst;
  Shared<ing> f;
  vec<int> order;
  KsumLookUpTable<ing> ksumtable(len);
  // Initialization sets all the pointers to nullptr inside.


  if (verbose) Rcpp::Rcout <<
    "Transform string input to bit streams and create flsss objects to mine..\n\n";


  // flsss decomposition will never use ksumtable. Pointers inside are all nullptr.
  bool ksumIsComputedInternallyAndWillNotBeExported = givenKsumTable.size() == 0;
  // If ksumIsComputedInternallyAndWillNotBeExported is true, then converting
  // strings to 64-bit integers, the container will be upsized to only
  // accommodate the largest subset sum rather than the superset sum.
  makeMflsssObjs<ing> (len, Xval, tsum, sizeNeeded, Ndecomposed, tp,
                       mflsssOBJvec, mitosisRst, f, tlimit, order, Mat, ksumtable,
                       ksumIsComputedInternallyAndWillNotBeExported);


  if (verbose)
  {
    Rcpp::Rcout << "The superset is compressed in a " << int(f.N) << " x " <<
      int(f.d) << " matrix of 64-bit integers." << "\n\n";


    Rcpp::Rcout << int(mflsssOBJvec.size()) <<
      " flsss objects are created and will be mined by " <<
        int(tp.maxCore) << " threads.\n\n";


    Rcpp::Rcout << int(mitosisRst.size()) <<
      " solutions are found during the decomposition.\n\n";
  }


  if (sizeNeeded - mitosisRst.size() > 0) // Mine flsss objects.
  {
    // Deal with ksumtable.
    processGivenKsumtableORrecompute<ing> (
        &mflsssOBJvec[0], mflsssOBJvec.size(), givenKsumTable, verbose,
        ksum, upscale, tp, ksumtable, f);


    if (verbose) Rcpp::Rcout << "Mining..\n\n";
    auto empFun = [](std::size_t t)->bool { return false; };


    tp.parFor(0, mflsssOBJvec.size(), [&](std::size_t i, std::size_t t)->bool
    {
      mflsssOBJvec[i].TTTstackRun();
      return f.totalSize >= f.sizeNeed or
        std::chrono::steady_clock::now() > f.endTime;
    }, 1, empFun, empFun);

  }


  int totalSol = mitosisRst.size();
  for(auto &x: mflsssOBJvec) totalSol += x.result.size();


  vec<vec<ing> > rst;
  rst.resize(totalSol);
  {
    int i = 0;
    for(int k = 0, kend = mitosisRst.size(); k < kend; ++k, ++i)
      rst[i].swap(mitosisRst[k]);
    for(int k = 0, kend = mflsssOBJvec.size(); k < kend; ++k)
    {
      for(int j = 0, jend = mflsssOBJvec[k].result.size(); j < jend; ++j, ++i)
        rst[i].swap(mflsssOBJvec[k].result[j]);
    }
  }


  List finalRst(rst.size());
  {
    for(int i = 0, iend = rst.size(); i < iend; ++i)
    {
      IntegerVector tmp(rst[i].size());
      for(int k = 0, kend = tmp.size(); k < kend; ++k)
        tmp[k] = order[rst[i][k]] + 1;
      finalRst[i] = tmp;
    }
  }


  return finalRst;
}




template<typename ing>
inline List mflsssArbObjRun(List X, int solutionNeed, double tlimit,
                            CharlieThreadPool &tp,
                            int ksum, int upscale, bool verbose)
{
  RawVector shared = X["shared"];
  Shared<ing> f;
  f.read(shared);
  IntegerVector order = X["order"];
  List trimat = X["trimat"];
  TriM Mat;
  Mat.read(trimat);
  f.sizeNeed = solutionNeed;
  f.totalSize = 0;
  f.endTime = std::chrono::steady_clock::now() +
    std::chrono::seconds(std::size_t(tlimit));
  f.M = Mat.mat;


  mflsssOBJ<ing> O;
  List xobj = X["obj"];
  O.read(xobj, &f);


  KsumLookUpTable<ing> ksumtable(f.subsetSize);
  Rcpp::List givenKsumTable = X["ksumtable"];


  processGivenKsumtableORrecompute<ing> (
      &O, 1, givenKsumTable, verbose, ksum, upscale, tp, ksumtable, f);
  tp.destroy();


  O.TTTstackRun();


  List rst(O.result.size());
  for (int i = 0, iend = O.result.size(); i < iend; ++i)
  {
    IntegerVector tmp(O.result[i].size());
    for (int k = 0, kend = tmp.size(); k < kend; ++k)
      tmp[k] = order[O.result[i][k]] + 1;
    rst[i] = tmp;
  }
  return rst;
}
















