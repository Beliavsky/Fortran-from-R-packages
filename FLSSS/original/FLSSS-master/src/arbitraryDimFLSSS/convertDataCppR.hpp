#pragma once
#include <Rcpp.h>


#ifndef vec
#define vec std::vector
#endif


inline std::uintptr_t ptrInt(void *x)
{
  return reinterpret_cast<std::uintptr_t> (x);
}


template<typename T> // pointer type is T.
inline T intPtr(std::uintptr_t x)
{
  return reinterpret_cast<T> (x);
}


template<typename T>
inline Rcpp::RawVector copy2rRaw(T &x)
{
  Rcpp::RawVector rst(sizeof(T));
  std::memcpy(&rst[0], &x, sizeof(T));
  return rst;
}


// template<typename T>
// inline Rcpp::RawVector copy2rRaw(T &&x)
// {
//   Rcpp::RawVector rst(sizeof(T));
//   std::memcpy(&rst[0], &x, sizeof(T));
//   return rst;
// }


template<typename T>
inline Rcpp::RawVector copyVec2rRaw(vec<T> &x)
{
  Rcpp::RawVector rst(x.size() * sizeof(T));
  std::memcpy(&rst[0], &x[0], rst.size());
  return rst;
}


template<typename T>
inline void copyRraw(T &x, Rcpp::RawVector v)
{
  std::memcpy((char*)(&x), (char*)(&v[0]), sizeof(T));
}


template<typename T>
inline T copyRraw(Rcpp::RawVector v)
{
  T x;
  std::memcpy((char*)(&x), (char*)(&v[0]), sizeof(T));
  return x;
}


template<typename T>
inline void copyRraw2vec(vec<T> &x, Rcpp::RawVector v)
{
  x.resize(v.size() / sizeof(T));
  std::copy(v.begin(), v.end(), (char*)(&x[0]));
}


// template<typename T>
// void zeroVecH(vec<T> &x)
// {
//   std::fill((char*)(&x[0]), (char*)(&x[0]) + sizeof(x), 0);
// }










