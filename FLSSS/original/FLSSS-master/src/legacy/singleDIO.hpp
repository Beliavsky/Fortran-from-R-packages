# pragma once
# include <Rcpp.h>
# include "macros.hpp"
# include "triMatOneD.hpp"
# include "PATclass.hpp"




template<typename valtype, typename indtype>
inline void importSKvec(
    void *buffer, std::size_t bufferSize, vec<PAT<valtype, indtype> > &SK,
    indtype *intCtnr,
    indtype *intCtnrPrior)
{
  SK.resize(bufferSize / sizeof(PAT<valtype, indtype>));
  std::memcpy(&SK[0], buffer, bufferSize);
  for(unsigned i = 0, iend = SK.size(); i < iend; ++i)
  {
    SK[i].LB = SK[i].LB - intCtnrPrior + intCtnr;
    SK[i].UB = SK[i].UB - intCtnrPrior + intCtnr;
    SK[i].UBleftResv = SK[i].UBleftResv - intCtnrPrior + intCtnr;
  }
}




template<typename valtype, typename indtype>
inline void importTriMat(void *buffer, std::size_t bufferSize,
                  triMoneD<valtype, indtype> &mat,
                  vec<unsigned char> &matContainer,
                  std::size_t L, std::size_t N)
{
  matContainer.resize(bufferSize);
  std::memcpy(&matContainer[0], buffer, bufferSize);
  mat.alloc(&matContainer[0], L, N);
}




template<typename valtype, typename indtype>
inline Rcpp::List exportAll(
    indtype len, indtype N, std::size_t SKbackOffset, valtype ME,
    vec<unsigned char> &matContainer, vec<indtype> &intCtnr,
    vec<PAT<valtype, indtype> > &SK, bool useBiSearchInFindBounds)
{
  Rcpp::IntegerVector len_N_SKback(4);
  len_N_SKback[0] = len;
  len_N_SKback[1] = N;
  len_N_SKback[2] = SKbackOffset;
  len_N_SKback[3] = useBiSearchInFindBounds;


  Rcpp::RawVector me(sizeof(valtype));
  std::memcpy(&me[0], &ME, sizeof(valtype));


  Rcpp::RawVector intCtnrAddress(sizeof(indtype*));
  indtype *tmpAddress = &intCtnr[0];
  std::memcpy(&intCtnrAddress[0], &tmpAddress, intCtnrAddress.size());


  Rcpp::RawVector matCt(matContainer.size());
  std::memcpy(&matCt[0], &matContainer[0], matCt.size());


  Rcpp::RawVector intCt(intCtnr.size() * sizeof(indtype));
  std::memcpy(&intCt[0], &intCtnr[0], intCt.size());


  Rcpp::RawVector SKct(SK.size() * sizeof(PAT<valtype, indtype>));
  std::memcpy(&SKct[0], &SK[0], SKct.size());


  Rcpp::IntegerVector intValSize(2);
  intValSize[0] = sizeof(indtype);
  intValSize[1] = sizeof(valtype);


  return Rcpp::List::create(
    Rcpp::Named("len_N_SKback") = len_N_SKback,
    Rcpp::Named("ME") = me,
    Rcpp::Named("intCtnrAddress") = intCtnrAddress,
    Rcpp::Named("matContainer") = matCt,
    Rcpp::Named("intCtnr") = intCt,
    Rcpp::Named("SK") = SKct,
    Rcpp::Named("types") = intValSize);
}




template<typename valtype, typename indtype>
inline void importAll(
    Rcpp::List memImage, indtype &len, indtype &N, valtype **&M, PAT<valtype, indtype> *&SKback, valtype &ME,
    vec<unsigned char> &matContainer, vec<indtype> &intCtnr,
    vec<PAT<valtype, indtype> > &SK, bool &useBiSearchInFindBounds)
{
  Rcpp::IntegerVector len_N_SKback = memImage[0];
  len = len_N_SKback[0];
  N = len_N_SKback[1];
  std::size_t SKbackOffset = len_N_SKback[2];
  useBiSearchInFindBounds = len_N_SKback[3];


  Rcpp::RawVector me = memImage[1];
  std::memcpy(&ME, &me[0], sizeof(valtype));


  Rcpp::RawVector memAddress = memImage[2];
  indtype *priorIntCnt;
  std::memcpy(&priorIntCnt, &memAddress[0], sizeof(indtype*));


  Rcpp::RawVector matCnt = memImage[3];
  triMoneD<valtype, indtype> mat;
  importTriMat<valtype, indtype> (&matCnt[0], matCnt.size(), mat, matContainer, len, N); // matContainer is resized inside function
  M = mat.mat;


  Rcpp::RawVector intCnt = memImage[4];
  intCtnr.resize(intCnt.size() / sizeof(indtype));
  std::memcpy(&intCtnr[0], &intCnt[0], intCnt.size());


  Rcpp::RawVector SKcnt = memImage[5];
  importSKvec<valtype, indtype> (&SKcnt[0], SKcnt.size(), SK, &intCtnr[0], priorIntCnt);
  SKback = &SK[0] + SKbackOffset;
}



