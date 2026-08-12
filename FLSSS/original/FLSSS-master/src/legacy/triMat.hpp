# pragma once
# include "mvalOperation.hpp"


template<typename valtype, typename indtype>
struct triM
{
  valtype ***mat; // mat is the vector of vector of vector


  std::size_t containerByteSize(std::size_t d, std::size_t N, std::size_t L)
  {
    std::size_t NofVal = (2 * N - L + 1) * L / 2;
    return NofVal * d * sizeof(valtype) + NofVal * sizeof(valtype*) + L * sizeof(valtype**);
  }


  inline void alloc(void *begin, indtype d, indtype len, indtype vlen)
  {
    std::size_t N = vlen, L = len;
    std::size_t NofMval = (2 * N - L + 1) * L / 2;
    valtype *contentEnd = (valtype*)begin + NofMval * d; // [begin, contentEnd) stores values


    valtype *mval = (valtype*)begin;
    valtype **mvalArray = (valtype**)contentEnd, **mvalArrayEnd = mvalArray + NofMval; // [mvalArray, mvalArrayEnd) stores pointers to mvals


    while(mvalArray < mvalArrayEnd)
    {
      *mvalArray = mval;
      ++mvalArray;
      mval += d;
    }
    mvalArray = (valtype **)contentEnd;
    valtype ***M = (valtype ***)mvalArrayEnd, ***Mend = M + L; // [M, Mend) stores pointers to the "columns" of M


    mat = M;
    indtype t = vlen;
    while(M < Mend)
    {
      *M = mvalArray;
      ++M;
      mvalArray += t;
      --t;
    }
  }


  inline valtype ** & operator [] (indtype i) {return mat[i];}


  inline void make(void *containerBegin, indtype len, Rcpp::List mv) // mv is the data frame of d (dimensions) columns
  {
    indtype d = mv.size();
    Rcpp::NumericVector tmp = mv[0];
    indtype vlen = tmp.size();
    alloc(containerBegin, d, len, vlen);


    // get the first col
    valtype **firstCol = mat[0];
    for(indtype k = 0, &kend = d; k < kend; ++k)
    {
      Rcpp::NumericVector v = mv[k];
      for(indtype i = 0, iend = vlen; i < iend; ++i)
      {
        firstCol[i][k] = v[i];
      }
    }


    valtype **v = mat[0];
    for(indtype i = 1, iend = len; i < iend; ++i)
    {
      valtype **lastCol = mat[i - 1], **currentCol = mat[i];
      for(indtype j = 0, jend = vlen - i; j < jend; ++j)
      {
        mvalPlus(currentCol[j], lastCol[j], v[i + j], d);
      }
    }
  }
};



