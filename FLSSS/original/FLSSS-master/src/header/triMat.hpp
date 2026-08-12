# pragma once
# include "mvalOperation.hpp"
# include "macros.hpp"




template<typename valtype, typename indtype>
struct triM
{
  valtype ***mat; // mat is the vector of vector of vector


  std::size_t containerWordSize(std::size_t d, std::size_t N, std::size_t L)
  {
    std::size_t NofVal = (2 * N - L + 1) * L / 2;
    std::size_t byteSizeAtLeast =
      NofVal * d * sizeof(valtype) + (NofVal + 1) *
      sizeof(valtype*) + (L + 1) * sizeof(valtype**);
    return byteSizeAtLeast / sizeof(word) + 1;
  }


  inline void alloc(void *begin, indtype d, indtype len, indtype vlen)
  {
    std::size_t N = vlen, L = len;
    std::size_t NofMval = (2 * N - L + 1) * L / 2;
    valtype *contentEnd = (valtype*)begin + NofMval * d; // [begin, contentEnd) stores values


    valtype *mval = (valtype*)begin;
    valtype **mvalArray = properAddress<valtype*>(contentEnd), **mvalArrayEnd = mvalArray + NofMval;
    // [mvalArray, mvalArrayEnd) stores pointers to mvals


    valtype **mvalArrayResv = mvalArray;
    while(mvalArray < mvalArrayEnd)
    {
      *mvalArray = mval;
      ++mvalArray;
      mval += d;
    }


    mvalArray = mvalArrayResv;
    valtype ***M = properAddress<valtype**>(mvalArrayEnd), ***Mend = M + L;


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


  inline valtype ** & operator [] (indtype i) { return mat[i]; }


  inline void make(void *containerBegin, indtype len, Rcpp::NumericMatrix mv)
    // mv is the matrix of d (dimensions) columns
  {
    indtype d = mv.ncol();
    indtype vlen = mv.nrow();
    alloc(containerBegin, d, len, vlen);


    // get the first col
    valtype **firstCol = mat[0];
    valtype *tmpv = (valtype*)&mv[0];
    for(indtype k = 0, &kend = d; k < kend; ++k)
    {
      for(indtype i = 0, iend = vlen; i < iend; ++i)
      {
        firstCol[i][k] = tmpv[i];
        // std::bitset<64> tb(firstCol[i][k]);
        // std::cout << tb << ", ";
      }
      // std::cout << "\n";
      tmpv += vlen;
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



