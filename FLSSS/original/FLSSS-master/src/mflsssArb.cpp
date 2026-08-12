// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]
#include "arbitraryDimFLSSS/mflsssDecomp.hpp"
using namespace Rcpp;


// [[Rcpp::export]]
List arbFLSSS(
    int len, StringMatrix V, StringVector target,
    SEXP givenKsumTable = R_NilValue,
    int solutionNeed = 1, int maxCore = 7, double tlimit = 60,
    int approxNinstance = 1000, int ksumK = 4, int ksumTableSizeScaler = 30,
    bool verbose = true)
{
  if (len <= 0)
  {
    stop("Subset size = 0 is currently unsupported. For mining subsets of \
arbitrary sizes, can the function on all subset sizes.\n");
  }


  maxCore = std::max<int> (1, maxCore);
  if (maxCore == 1) approxNinstance = 1;


  CharlieThreadPool tp(maxCore);


  if (TYPEOF(givenKsumTable) == NILSXP) givenKsumTable = List::create();
  int maxint = std::max(V.ncol(), V.nrow());
  List rst;
  if (maxint < 127) rst = mflsssArb<int8_t> (
    len, V, target, solutionNeed, tp, approxNinstance, tlimit,
    givenKsumTable, ksumK, ksumTableSizeScaler, verbose);
  else if (maxint < 32767) rst = mflsssArb<int16_t> (
    len, V, target, solutionNeed, tp, approxNinstance, tlimit,
    givenKsumTable, ksumK, ksumTableSizeScaler, verbose);
  else rst = mflsssArb<int32_t> (
    len, V, target, solutionNeed, tp, approxNinstance, tlimit,
    givenKsumTable, ksumK, ksumTableSizeScaler, verbose);
  return rst;
}


// [[Rcpp::export]]
List decomposeArbFLSSS(
    int len, StringMatrix V, StringVector target,
    int approxNinstance = 1000, int maxCore = 7,
    SEXP ksumTable = R_NilValue,
    int ksumK = 4, int ksumTableSizeScaler = 30, bool verbose = true)
{
  if (len <= 0)
  {
    stop("Subset size = 0 is currently unsupported. For mining subsets of \
           arbitrary sizes, can the function on all subset sizes.\n");
  }


  int maxint = std::max(V.nrow(), V.ncol());
  List rst;
  CharlieThreadPool tp(maxCore);
  if (maxint < 127) rst = decomposeArbFLSSSexport<int8_t> (
    len, V, target, approxNinstance, tp, ksumTable, ksumK,
    ksumTableSizeScaler, verbose);
  else if (maxint < 32767) rst = decomposeArbFLSSSexport<int16_t> (
    len, V, target, approxNinstance, tp, ksumTable, ksumK,
    ksumTableSizeScaler, verbose);
  else rst = decomposeArbFLSSSexport<int32_t> (
    len, V, target, approxNinstance, tp, ksumTable, ksumK,
    ksumTableSizeScaler, verbose);
  return rst;
}


// [[Rcpp::export]]
List arbFLSSSobjRun(List X, int solutionNeed = 1,
                    double tlimit = 60, int maxCore = 7, int ksumK = 0,
                    int ksumTableSizeScaler = 30, bool verbose = true)
{
  if (X.size() == 0) return List::create();
  int indtype = X["indtype"];
  List rst;
  CharlieThreadPool tp(maxCore);
  if (indtype == 1)
    rst = mflsssArbObjRun<int8_t> (
      X, solutionNeed, tlimit, tp, ksumK, ksumTableSizeScaler, verbose);
  else if (indtype == 2)
    rst = mflsssArbObjRun<int16_t> (
      X, solutionNeed, tlimit, tp, ksumK, ksumTableSizeScaler, verbose);
  else if (indtype == 4)
    rst = mflsssArbObjRun<int32_t> (
      X, solutionNeed, tlimit, tp, ksumK, ksumTableSizeScaler, verbose);
  else
    rst = mflsssArbObjRun<int64_t> (
      X, solutionNeed, tlimit, tp, ksumK, ksumTableSizeScaler, verbose);
  return rst;
}


// lb and ub's sizes must be no less than ksum.
List ksumHashNoTarget(int ksumK, uint64_t *x, int d, int N,
                      int *lb, int *ub, int upscale,
                      bool verbose, CharlieThreadPool &tp)
{
  Ksum KS;
  Rcpp::List rst(std::max<int> (ksumK - 3 + 1, 0));
  Rcpp::StringVector enames(rst.size());
  Rcpp::StringVector eleName(2);
  eleName[0] = "prime"; eleName[1] = "table";


  for (int k = 3; k <= ksumK; ++k)
  {
    if (verbose)
    {
      Rcout << "Index lower bounds for " << k << "-sum: ";
      for (int u = 0; u < k; ++u) Rcout << lb[u] << ", ";
      Rcout << "\n";
      Rcout << "Index Upper bounds for " << k << "-sum: ";
      for (int u = 0; u < k; ++u) Rcout << ub[ksumK - k + u] << ", ";
      Rcout << "\n";
    }

    KS.reset(x, N, d, k, lb, ub + (ksumK - k), upscale, tp);
    RawVector tmp = KS(verbose);
    RawVector tmpprime = copy2rRaw<uint64_t> (KS.modPrime);
    auto elelist = Rcpp::List::create(tmpprime, tmp);
    elelist.names() = eleName;
    rst[k - 3] = elelist;
    enames[k - 3] = std::to_string(k) + "-sum";
  }
  rst.names() = enames;


  return rst;
}


// [[Rcpp::export]]
List ksumHash(int ksumK, StringMatrix V, int ksumTableSizeScaler = 30,
              SEXP target = R_NilValue, int len = 0,
              int approxNinstance = 1000, bool verbose = true, int maxCore = 7)
{
  bool ksumIsComputedInternallyAndWillNotBeExported = false;
  ksumK = std::min<int> (std::max<int> (3, ksumK), V.nrow());


  CharlieThreadPool tp(maxCore);


  if (TYPEOF(target) == NILSXP)
  {
    auto targetNew = StringVector(V.ncol(), "0");
    if (verbose) Rcpp::Rcout <<
      "Transform string input to bit streams..\n\n";
    vec<uint64_t> Xbit, targetSS, largestSS;
    vec<int> colOrder, order;
    int totalBits = 0;
    stringMatTo64bitIntMat(
      V, targetNew, 1, tp, Xbit, targetSS, largestSS,
      order, colOrder, totalBits,
      ksumIsComputedInternallyAndWillNotBeExported);
    int d = targetSS.size(), N = Xbit.size() / d;


    vec<int> lub(ksumK * 2);
    for (int i = 0; i < ksumK; ++i)
    {
      lub[i] = i;
      lub[i + ksumK] = V.nrow() - ksumK + i;
    }


    // for (int i = 0; i < ksumK; ++i) Rcpp::Rcout << lub[i] << ", ";
    // Rcpp::Rcout << "\n";
    // for (int i = 0; i < ksumK; ++i) Rcpp::Rcout << lub[i + ksumK] << ", ";
    // Rcpp::Rcout << "\n";


    return ksumHashNoTarget(ksumK, &Xbit[0], d, N, &lub[0], &lub[ksumK],
                            ksumTableSizeScaler, verbose, tp);
  }


  ksumK = std::min<int> (std::max<int> (3, ksumK), len);
  vec<mflsssOBJ<int> > mflsssOBJvec;
  TriM Mat;
  vec<vec<int> > mitosisRst;
  vec<int> order;
  KsumLookUpTable<int> ksumtableDummy(len);
  if (verbose) Rcpp::Rcout <<
    "Transform string input to bit streams and create flsss objects..\n\n";
  Shared<int> fdummy;
  if (len <= 0) stop("Subset size = 0 is currently unsupported.\n");


  makeMflsssObjs<int> (len, V, target, 1e9, approxNinstance, tp,
                       mflsssOBJvec, mitosisRst, fdummy, 1e9, order, Mat,
                       ksumtableDummy,
                       ksumIsComputedInternallyAndWillNotBeExported);


  if (verbose)
  {
    Rcpp::Rcout << int(mflsssOBJvec.size()) << " flsss objects are created.\n\n";
    Rcpp::Rcout << int(mitosisRst.size()) <<
      " solutions are found during the decomposition.\n\n";
  }


  List ksumTableRst(0);
  processGivenKsumtableORrecompute<int> (
      &mflsssOBJvec[0], mflsssOBJvec.size(), ksumTableRst, verbose,
      ksumK, ksumTableSizeScaler, tp, ksumtableDummy, fdummy);


  return ksumTableRst;
}










/*
 // [[Rcpp::export]]
 NumericVector ifAddrSame(List X)
 {
 NumericVector rst(X.size());
 for (int i = 0, iend = rst.size(); i < iend; ++i)
 {
 List tmplist  = X[i];
 IntegerVector tmp = tmplist["order"];
 rst[i] = (double)(std::uintptr_t(&tmp[0]));
 }
 return rst;
 }
 */


































