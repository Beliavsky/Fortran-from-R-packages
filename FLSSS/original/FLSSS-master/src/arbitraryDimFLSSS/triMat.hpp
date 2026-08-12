#pragma once
#include "arithmetic.hpp"
#include <Rcpp.h>
#include "convertDataCppR.hpp"


// DO NOT DELETE!
struct TriM
{
  uint64_t ***mat; // mat is the vector of vector of vector
  Rcpp::RawVector content;


  void reset(uint64_t *x, uint64_t d, uint64_t N, uint64_t L)
  {
    uint64_t Nval = (2 * N - L + 1) * L / 2;
    uint64_t &NdoublePtr = L, &NsinglePtr = Nval;
    uint64_t bytes = Nval * sizeof(uint64_t) * d +
      NdoublePtr * sizeof(uint64_t**) + sizeof(uint64_t**) * Nval + 64;
    content = Rcpp::RawVector(bytes);
    mat = (uint64_t***) N8BAA(&content[0]);
    uint64_t **singlePtrBegin = (uint64_t**) N8BAA(mat + NdoublePtr);
    uint64_t *valBegin = (uint64_t*) N8BAA(singlePtrBegin + NsinglePtr);


    // Fill single pointers.
    for(uint64_t i = 0; i < Nval;
      singlePtrBegin[i] = valBegin + i * d, ++i);


    // Fill double pointers;
    for(uint64_t i = 0, **x = singlePtrBegin; i < L;
      mat[i] = x, x += N - i, ++i);


    std::memcpy(mat[0][0], x, sizeof(uint64_t) * N * d);
    uint64_t **v = mat[0];
    for(int64_t i = 1, iend = L; i < iend; ++i)
    {
      uint64_t **lastCol = mat[i - 1], **currentCol = mat[i];
      for(int64_t j = 0, jend = N - i; j < jend; ++j)
        mvalPlus(currentCol[j], lastCol[j], v[i + j], d);
    }


    // // Better not DELETE IT BABY.
    // std::cout << "tri mat = " << std::endl;
    // for(int i = 0, iend = L; i < iend; ++i)
    // {
    //   for(int j = 0, jend = N - i; j < jend; ++j)
    //   {
    //     for(int k = 0, kend = d; k < kend; ++k)
    //       std::cout << mat[i][j][k] << ", ";
    //     std::cout << std::endl;
    //   }
    //   std::cout << std::endl;
    // }
  }


  TriM(){}
  TriM(uint64_t *x, uint64_t d, uint64_t N, uint64_t L) { reset(x, d, N, L); }


  Rcpp::List save()
  {
    auto tmpPtrInt = ptrInt(mat);
    Rcpp::RawVector contentAddr = copy2rRaw<std::uintptr_t> (tmpPtrInt);
    tmpPtrInt = ptrInt(&mat[0][0][0]);
    Rcpp::RawVector firstValAddr = copy2rRaw<std::uintptr_t> (tmpPtrInt);
    return Rcpp::List::create(
      Rcpp::Named("contentAddr") = contentAddr,
      Rcpp::Named("firstValAddr") = firstValAddr,
      Rcpp::Named("content") = content);
  }


  void read(Rcpp::List trimat)
  {
    Rcpp::RawVector contentAddrRaw = trimat["contentAddr"];
    Rcpp::RawVector firstValAddrRaw = trimat["firstValAddr"];
    content = trimat["content"];
    mat = (uint64_t***) N8BAA(&content[0]);
    std::uintptr_t contentAddr, firstValAddr;
    copyRraw<std::uintptr_t> (contentAddr, contentAddrRaw);
    copyRraw<std::uintptr_t> (firstValAddr, firstValAddrRaw);
    std::uintptr_t offset = ptrInt(mat) - contentAddr;
    std::uintptr_t Nptr = (firstValAddr - contentAddr) / sizeof(std::uintptr_t);
    std::uintptr_t *p = (std::uintptr_t*)mat;
    for (std::uintptr_t i = 0; i < Nptr; ++i) p[i] += offset;
  }


};









