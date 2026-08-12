# pragma once
# include "oneDoperation.hpp"
# include <fstream>


namespace legacy
{
template<typename valtype, typename indtype>
struct triMoneD
{
  valtype **mat; // mat is the vector of vector of vector


  std::size_t containerByteSize(std::size_t N, std::size_t L)
  {
    std::size_t NofVal = (2 * N - L + 1) * L / 2;
    return NofVal * sizeof(valtype) + L * sizeof(valtype*);
  }


  inline void alloc(void *begin, indtype len, indtype vlen)
  {
    std::size_t N = vlen, L = len;
    std::size_t NofMval = (2 * N - L + 1) * L / 2;
    valtype *contentEnd = (valtype*)begin + NofMval; // [begin, contentEnd) stores values


    mat = (valtype**)contentEnd;
    mat[0] = (valtype*)begin;
    for(indtype t = 1; t < len; ++t)
    {
      mat[t] = mat[t - 1] + (vlen - t + 1);
    }
  }


  void print(int N, int L)
  {
    std::ofstream outfile("triMat.csv");
    for(int i = 0; i < L; ++i)
    {
      for(int j = 0, jend = N - i; j < jend; ++j)
      {
        outfile << mat[i][j] << ",";
      }
      outfile << "\n";
    }
  }


  inline void make(void *containerBegin, indtype len, Rcpp::NumericVector v) // mv is the data frame of d (dimensions) columns
  {
    indtype vlen = v.size();
    alloc(containerBegin, len, vlen);


    // get the first col
    // std::copy(v.begin(), v.end(), mat[0]);
    std::memcpy(mat[0], &*v.begin(), sizeof(valtype) * v.size());


    valtype *V = mat[0];
    for(indtype i = 1, iend = len; i < iend; ++i)
    {
      valtype *lastCol = mat[i - 1], *currentCol = mat[i];
      for(indtype j = 0, jend = vlen - i; j < jend; ++j)
      {
        currentCol[j] = lastCol[j] + V[i + j];
      }
    }
  }


};
}


