# pragma once
# include <Rcpp.h>
# include "macros.hpp"
# include "mPATclass.hpp"
# include "triMat.hpp"


template<typename valtype, typename indtype>
inline void importSKvec(
    void *buffer, std::size_t bufferSize, vec<mPAT<valtype, indtype> > &SK,
    indtype *intCtnr, valtype *valCtnr,
    indtype *intCtnrPrior, valtype *valCtnrPrior)
{
  SK.resize(bufferSize / sizeof(mPAT<valtype, indtype>));
  std::memcpy(&SK[0], buffer, bufferSize);
  for(unsigned i = 0, iend = SK.size(); i < iend; ++i)
  {
    SK[i].LB = SK[i].LB - intCtnrPrior + intCtnr;
    SK[i].UB = SK[i].UB - intCtnrPrior + intCtnr;
    SK[i].Bresv = SK[i].Bresv - intCtnrPrior + intCtnr;
    SK[i].MIN = SK[i].MIN - valCtnrPrior + valCtnr;
    SK[i].MAX = SK[i].MAX - valCtnrPrior + valCtnr;


    SK[i].sumLB = SK[i].sumLB - valCtnrPrior + valCtnr;
    SK[i].sumUB = SK[i].sumUB - valCtnrPrior + valCtnr;
  }
}




template<typename valtype, typename indtype>
inline void importTriMat(void *buffer, std::size_t bufferSize,
                  triM<valtype, indtype> &mat,
                  vec<unsigned char> &matContainer,
                  std::size_t d, std::size_t L, std::size_t N)
{
  matContainer.resize(bufferSize);
  std::memcpy(&matContainer[0], buffer, bufferSize);
  mat.alloc(&matContainer[0], d, L, N);
}




// Single threaded version for mFLSSScomo
template<typename valtype, typename indtype>
inline Rcpp::List exportAll(
    indtype d, indtype dlst, indtype dl, indtype dust, indtype du,
    indtype len, indtype N, std::size_t SKbackOffset,
    vec<unsigned char> &matContainer, vec<indtype> &intCtnr, vec<valtype> &valCtnr,
    vec<mPAT<valtype, indtype> > &SK, bool useBiSearchInFindBounds)
{
  Rcpp::IntegerVector d_len_N_SKback(9);
  d_len_N_SKback[0] = d;
  d_len_N_SKback[1] = dlst;
  d_len_N_SKback[2] = dl;
  d_len_N_SKback[3] = dust;
  d_len_N_SKback[4] = du;
  d_len_N_SKback[5] = len;
  d_len_N_SKback[6] = N;
  d_len_N_SKback[7] = SKbackOffset;
  d_len_N_SKback[8] = useBiSearchInFindBounds;


  Rcpp::RawVector intCtnrValCtnrHeadAddress(sizeof(indtype *) + sizeof(valtype *));
  void *pc[2];
  pc[0] = &intCtnr[0];
  pc[1] = &valCtnr[0];
  std::memcpy(&intCtnrValCtnrHeadAddress[0], pc, intCtnrValCtnrHeadAddress.size());


  Rcpp::RawVector matCt(matContainer.size());
  std::memcpy(&matCt[0], &matContainer[0], matCt.size());


  Rcpp::RawVector intCt(intCtnr.size() * sizeof(indtype));
  std::memcpy(&intCt[0], &intCtnr[0], intCt.size());


  Rcpp::RawVector valCt(valCtnr.size() * sizeof(valtype));
  std::memcpy(&valCt[0], &valCtnr[0], valCt.size());


  Rcpp::RawVector SKct(SK.size() * sizeof(mPAT<valtype, indtype>));
  std::memcpy(&SKct[0], &SK[0], SKct.size());


  Rcpp::IntegerVector intValSize(2);
  intValSize[0] = sizeof(indtype);
  intValSize[1] = sizeof(valtype);


  return Rcpp::List::create(
    Rcpp::Named("d_len_N_SKback") = d_len_N_SKback,
    // Rcpp::Named("ME") = me,
    Rcpp::Named("intCtnrValCtnrHeadAddress") = intCtnrValCtnrHeadAddress,
    Rcpp::Named("matContainer") = matCt,
    Rcpp::Named("intCtnr") = intCt,
    Rcpp::Named("valCtnr") = valCt,
    Rcpp::Named("SK") = SKct,
    Rcpp::Named("types") = intValSize);
}




// Single threaded version for mFLSSScomo
template<typename valtype, typename indtype>
inline void importAll(Rcpp::List memImage,
    indtype &d, indtype &dlst, indtype &dl, indtype &dust, indtype &du,
    indtype &len, indtype &N, valtype ***&M, mPAT<valtype, indtype> *&SKback,
    vec<unsigned char> &matContainer, vec<indtype> &intCtnr,
    vec<valtype> &valCtnr, vec<mPAT<valtype, indtype> > &SK, bool &useBiSearchInFindBounds)
{
  Rcpp::IntegerVector d_len_N_SKback = memImage[0];


  d = d_len_N_SKback[0];
  dlst = d_len_N_SKback[1];
  dl = d_len_N_SKback[2];
  dust = d_len_N_SKback[3];
  du = d_len_N_SKback[4];
  len = d_len_N_SKback[5];
  N = d_len_N_SKback[6];
  std::size_t SKbackOffset = d_len_N_SKback[7];
  useBiSearchInFindBounds = d_len_N_SKback[8];


  Rcpp::RawVector memAddress = memImage[1];
  indtype *priorIntCnt;
  valtype *priorValCnt;
  std::memcpy(&priorIntCnt, &memAddress[0], sizeof(indtype*));
  std::memcpy(&priorValCnt, &memAddress[0] + sizeof(indtype*), sizeof(valtype*));


  Rcpp::RawVector matCnt = memImage[2];
  triM<valtype, indtype> mat;
  importTriMat<valtype, indtype> (
      &matCnt[0], matCnt.size(), mat, matContainer, d, len, N); // matContainer is resized inside function
  M = mat.mat;


  Rcpp::RawVector intCnt = memImage[3];
  intCtnr.resize(intCnt.size() / sizeof(indtype));
  std::memcpy(&intCtnr[0], &intCnt[0], intCnt.size());


  Rcpp::RawVector vanCnt = memImage[4];
  valCtnr.resize(vanCnt.size() / sizeof(valtype));
  std::memcpy(&valCtnr[0], &vanCnt[0], vanCnt.size());


  Rcpp::RawVector SKcnt = memImage[5];
  importSKvec<valtype, indtype> (&SKcnt[0], SKcnt.size(), SK, &intCtnr[0], &valCtnr[0], priorIntCnt, priorValCnt);
  SKback = &SK[0] + SKbackOffset;
}








// ________________________________________________________________________________________________
// Multi threaded version for mFLSSS
/*
template<typename valtype, typename indtype>
inline Rcpp::List exportAll(
    indtype d, indtype len, indtype N, int targetStartIndex,
    vec<std::size_t> &SKbackOffset, valtype *ME,
    vec<unsigned char> &matContainer,
    vec<vec<indtype> > &intCtnrGroup, vec<vec<valtype> > &valCtnrGroup,
    vec<vec<mPAT<valtype, indtype> > > &SKgroup,
    Rcpp::NumericVector keyTarget, Rcpp::NumericVector scaleFactor,
    Rcpp::NumericVector originalTarget,
    bool useBiSearchInFindBounds)
{
  Rcpp::IntegerVector d_len_N(5);
  d_len_N[0] = d;
  d_len_N[1] = len;
  d_len_N[2] = N;
  d_len_N[3] = useBiSearchInFindBounds;
  d_len_N[4] = targetStartIndex;


  Rcpp::IntegerVector SKbackOffsetGroup(SKbackOffset.begin(), SKbackOffset.end());


  Rcpp::NumericVector me(ME, ME + d);


  Rcpp::RawVector intValCntrHeadAddressGroup(
      (sizeof(indtype *) + sizeof(valtype *)) * intCtnrGroup.size());
  {
    char *st = (char *)&intValCntrHeadAddressGroup[0];
    unsigned siz = sizeof(indtype *) + sizeof(valtype *);
    for(int i = 0, iend = intCtnrGroup.size(); i < iend; ++i)
    {
      void *pc[2];
      pc[0] = &intCtnrGroup[i][0];
      pc[1] = &valCtnrGroup[i][0];
      std::memcpy(st, pc, siz);
      st += siz;
    }
  }


  Rcpp::RawVector matCt(matContainer.size());
  std::memcpy(&matCt[0], &matContainer[0], matCt.size());


  Rcpp::List intGroup(intCtnrGroup.size());
  {
    for(int i = 0, iend = intCtnrGroup.size(); i < iend; ++i)
    {
      Rcpp::RawVector tmp(intCtnrGroup[i].size() * sizeof(indtype));
      std::memcpy(&tmp[0], &intCtnrGroup[i][0], tmp.size());
      intGroup[i] = tmp;
    }
  }


  Rcpp::List valGroup(valCtnrGroup.size());
  {
    for(int i = 0, iend = valGroup.size(); i < iend; ++i)
    {
      Rcpp::RawVector tmp(valCtnrGroup[i].size() * sizeof(valtype));
      std::memcpy(&tmp[0], &valCtnrGroup[i][0], tmp.size());
      valGroup[i] = tmp;
    }
  }


  Rcpp::List stackGroup(SKgroup.size());
  {
    for(int i = 0, iend = stackGroup.size(); i < iend; ++i)
    {
      Rcpp::RawVector tmp(SKgroup[i].size() * sizeof(mPAT<valtype, indtype>));
      std::memcpy(&tmp[0], &SKgroup[i][0], tmp.size());
      stackGroup[i] = tmp;
    }
  }


  Rcpp::IntegerVector intValSize(2);
  intValSize[0] = sizeof(indtype);
  intValSize[1] = sizeof(valtype);


  return Rcpp::List::create(
    Rcpp::Named("d_len_N") = d_len_N,
    Rcpp::Named("SKbackOffsetGroup") = SKbackOffsetGroup,
    Rcpp::Named("ME") = me,
    Rcpp::Named("intValCntrHeadAddressGroup") = intValCntrHeadAddressGroup,
    Rcpp::Named("matContainer") = matCt,
    Rcpp::Named("intGroup") = intGroup,
    Rcpp::Named("valGroup") = valGroup,
    Rcpp::Named("stackGroup") = stackGroup,
    Rcpp::Named("originalTarget") = originalTarget,
    Rcpp::Named("keyTarget") = keyTarget,
    Rcpp::Named("scaleFactor") = scaleFactor,
    Rcpp::Named("types") = intValSize);
}




// I really don't have the energy for crafting a more comprehensive export function for multithreaded version now.
template<typename valtype, typename indtype>
inline Rcpp::List exportAll(
    indtype d, indtype len, indtype N, int firstRunningTask,
    valtype *ME, vec<unsigned char> &matContainer,
    Rcpp::NumericVector keyTarget, Rcpp::NumericVector scaleFactor,
    Rcpp::NumericVector originalTarget, Rcpp::List priorResult,
    bool useBiSearchInFindBounds)
{
  Rcpp::IntegerVector d_len_N(4);
  d_len_N[0] = d;
  d_len_N[1] = len;
  d_len_N[2] = N;
  d_len_N[3] = useBiSearchInFindBounds;


  Rcpp::NumericVector me(ME, ME + d);


  Rcpp::RawVector matCt(matContainer.size());
  std::memcpy(&matCt[0], &matContainer[0], matCt.size());


  Rcpp::IntegerVector intValSize(2);
  intValSize[0] = sizeof(indtype);
  intValSize[1] = sizeof(valtype);


  Rcpp::NumericVector remainKeyTarget(keyTarget.begin() + firstRunningTask, keyTarget.end());


  return Rcpp::List::create(
    Rcpp::Named("d_len_N") = d_len_N,
    Rcpp::Named("ME") = me,
    Rcpp::Named("matContainer") = matCt,
    Rcpp::Named("originalTarget") = originalTarget,
    Rcpp::Named("keyTarget") = remainKeyTarget,
    Rcpp::Named("scaleFactor") = scaleFactor,
    Rcpp::Named("types") = intValSize);
}




// Multi threaded version
template<typename valtype, typename indtype>
inline void importAll(
    Rcpp::List memImage,
    indtype &d, indtype &len, indtype &N, int &targetStartIndex, valtype ***&M,
    vec<valtype> &me, vec<unsigned char> &matContainer,
    Rcpp::NumericVector keyTarget, Rcpp::NumericVector scaleFactor,
    Rcpp::NumericVector originalTarget, bool &useBiSearchInFindBounds)
{
  Rcpp::IntegerVector d_len_N = memImage[0];
  d = d_len_N[0];
  len = d_len_N[1];
  N = d_len_N[2];
  useBiSearchInFindBounds = d_len_N[3];
  // targetStartIndex = d_len_N[4];


  Rcpp::NumericVector mer = memImage[1];
  me.assign(mer.begin(), mer.end());


  triM<valtype, indtype> mat;
  Rcpp::RawVector matCt = memImage[2];
  importTriMat<valtype, indtype> (&matCt[0], matCt.size(), mat, matContainer, d, len, N);
  M = mat.mat;


  originalTarget = memImage[3];
  keyTarget = memImage[4];
  scaleFactor = memImage[5];
}
*/



