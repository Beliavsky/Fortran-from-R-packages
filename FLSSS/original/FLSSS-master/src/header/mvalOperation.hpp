# pragma once
// # include <Rcpp.h>
# include "macros.hpp"
# include <bitset>
# include <fstream>
// using namespace Rcpp;


template<typename T>
inline T* properAddress(void *current)
{
  std::size_t adr = (std::size_t)current;
  std::size_t sizeofT = sizeof(T);
  std::size_t q = adr / sizeofT, r = adr % sizeofT;
  if(r == 0) return (T*)current;
  return (T*)(q * sizeofT) + 1;
}


template<typename T>
inline void subset(vec<T> &X, bool *B, unsigned worthItems)
{
  vec<T> Y(worthItems);
  unsigned j = 0;
  for(unsigned i = 0, iend = X.size(); i < iend; ++i)
  {
    if(B[i])
    {
      std::swap(X[i], Y[j]);
      ++j;
    }
  }
  std::swap(X, Y);
}


template<typename valtype, typename indtype>
inline void iterSum(valtype *rst, valtype **V, indtype *ele, indtype eleSize, indtype d)
{
  for(indtype i = 0; i < eleSize; ++i)
  {
    valtype* val = V[ele[i]];
    for(indtype k = 0; k < d; ++k)
    {
      rst[k] += val[k];
    }
  }
}


template<typename valtype, typename indtype>
inline void mvalPlus(valtype *tbegin, valtype *left, valtype *right, indtype d)
{
  for(valtype*&i = tbegin, *iend = i + d; i < iend; ++i, ++left, ++right)
  {
    *i = *left + *right;
  }
}


template<typename valtype, typename indtype>
inline void mvalMinus(valtype *tbegin, valtype *left, valtype *right, indtype d)
{
  for(valtype*&i = tbegin, *iend = i + d; i < iend; ++i, ++left, ++right)
  {
    *i = *left - *right;
  }
}


template<typename valtype, typename indtype>
inline void mvalPlusMinus(valtype *tbegin, valtype *left, valtype *mid, valtype *right, indtype d)
{
  for(valtype *&i = tbegin, *iend = i + d; i < iend; ++i, ++left, ++mid, ++right)
  {
    *i = *left + *mid - *right;
  }
}


template<typename valtype, typename indtype>
inline void mvalMinusPlus(valtype *tbegin, valtype *left, valtype *mid, valtype *right, indtype d)
{
  for(valtype*&i = tbegin, *iend = i + d; i < iend; ++i, ++left, ++mid, ++right)
  {
    *i = *left - *mid + *right;
  }
}




template<typename valtype, typename indtype, bool mk>
inline bool notAllLessEqual(valtype *x, valtype *y, indtype &i, indtype d, INT *mask)
{
  if(!mk)
  {
    for(; i < d; ++i)
    {
      if(x[i] > y[i]) return 1;
    }
  }
  else
  {
    for(; i < d; ++i)
    {
      // if: (y[i] - x[i]) & mask[i] == 0, then every one in y >= those in x
      // if: (y[i] - x[i]) & mask[i], then someone in x is greater than in y
      // {
      //   std::ofstream of("debug.csv", std::ios::app);
      //   std::bitset<64> tmp(x[i]);
      //   of << "\nx[i] = " << tmp << "\n";
      //   std::bitset<64> tmp2(y[i]);
      //   of << "(INT(y[i] - x[i]) & mask[i]) != 0 = " <<
      //     int((INT(y[i] - x[i]) & mask[i]) != 0) << "\n";
      //   of.close();
      // }
      if((INT(y[i] - x[i]) & mask[i]) != 0) return 1;
    }
  }
  // valtype would only be double, and it has the same bit size as INT
  return 0;
}


template<typename valtype, typename indtype, bool mk>
inline bool notAllLessEqual(valtype *x, valtype *y, indtype d, INT *mask)
{
  if(!mk)
  {
    for(indtype i = 0; i < d; ++i)
    {
      if(x[i] > y[i]) return 1;
    }
  }
  else
  {
    for(indtype i = 0; i < d; ++i)
    {
      // {
      //   std::ofstream of("debug.csv", std::ios::app);
      //   std::bitset<64> tmp(x[i]);
      //   of << "\nx[i] = " << tmp << "\n";
      //   std::bitset<64> tmp2(y[i]);
      //   of << "\ny[i] = " << tmp2 << "\n";
      //   of << "(INT(y[i] - x[i]) & mask[i]) != 0 = " <<
      //     int((INT(y[i] - x[i]) & mask[i]) != 0) << "\n";
      //   of.close();
      // }
      if((INT(y[i] - x[i]) & mask[i]) != 0) return 1;
    }
  }
  return 0;
}


template<typename valtype, typename indtype, bool mk>
inline bool notAllGreaterEqual(valtype *x, valtype *y, indtype &i, indtype d, INT *mask)
{
  if(!mk)
  {
    for(; i < d; ++i)
    {
      if(x[i] < y[i]) return 1;
    }
  }
  else
  {
    for(; i < d; ++i)
    {
      // {
      //   std::ofstream of("debug.csv", std::ios::app);
      //   std::bitset<64> tmp(x[i]);
      //   of << "\nx[i] = " << tmp << "\n";
      //   std::bitset<64> tmp2(y[i]);
      //   of << "\ny[i] = " << tmp2 << "\n";
      //   of << "(INT(x[i] - y[i]) & mask[i]) != 0 = " <<
      //     int((INT(x[i] - y[i]) & mask[i]) != 0) << "\n";
      //   of.close();
      // }
      if((INT(x[i] - y[i]) & mask[i]) != 0) return 1;
    }
  }
  return 0;
}


template<typename valtype, typename indtype, bool mk>
inline bool notAllGreaterEqual(valtype *x, valtype *y, indtype d, INT *mask)
{
  if(!mk)
  {
    for(indtype i = 0; i < d; ++i)
    {
      if(x[i] < y[i]) return 1;
    }
  }
  else
  {
    for(indtype i = 0; i < d; ++i)
    {
      // {
      //   std::ofstream of("debug.csv", std::ios::app);
      //   std::bitset<64> tmp(x[i]);
      //   of << "\nx[i] = " << tmp << "\n";
      //   std::bitset<64> tmp2(y[i]);
      //   of << "\ny[i] = " << tmp2 << "\n";
      //   of << "(INT(x[i] - y[i]) & mask[i]) != 0 = " <<
      //     int((INT(x[i] - y[i]) & mask[i]) != 0) << "\n";
      //   of.close();
      // }
      if((INT(x[i] - y[i]) & mask[i]) != 0) return 1;
    }
  }
  return 0;
}


// template<typename valtype, typename indtype>
// inline bool allEqual(valtype *x, valtype *y, indtype &i, indtype d)
// {
//   for(; i < d; ++i)
//   {
//     if(std::abs(x[i] - y[i]) > eps) return 0;
//   }
//   return 1;
// }
//
//
// template<typename valtype, typename indtype>
// inline bool allEqual(valtype *x, valtype *y, indtype d)
// {
//   for(indtype i = 0; i < d; ++i)
//   {
//     if(std::abs(x[i] - y[i]) > eps) return 0;
//   }
//   return 1;
// }








// ————————————————————————————————————————————————————————————————————————————————————————————————
template<typename valtype, typename indtype, bool mk> // this function is for experimenting
// inline valtype **mvalLowerBoundBiMan(valtype **begin, valtype **end, valtype *val, indtype d)
inline valtype **mvalLowerBoundBiMan(valtype **begin, valtype **end, valtype *val,
                                     indtype compst, indtype compsiz, INT *mask)
{
  // if(!notAllGreaterEqual(*begin, val, d)) return begin;
  if(!notAllGreaterEqual<valtype, indtype, mk> (
      *begin + compst, val + compst, compsiz, mask)) return begin;
  --end;
  while(true)
  {
    valtype **mid = begin + int((end - begin) / 2);
    // if(notAllGreaterEqual(*mid, val, d))
    if(notAllGreaterEqual<valtype, indtype, mk> (
        *mid + compst, val + compst, compsiz, mask))
    {
      if(mid == begin)
      {
        return end;
      }
      begin = mid;
    }
    // to here, it means mid is all greater equal than val
    // else if(!notAllGreaterEqual(*(mid - 1), val, d))
    else if(!notAllGreaterEqual<valtype, indtype, mk> (
        *(mid - 1) + compst, val + compst, compst, mask))
    // and at the same time even mid-1 is all greater equal than val
    {
      end = mid; // so end is safe to set as mid
    }
    // to here, it means mid is the first that all greater equal than val from low to high
    else
    {
      return mid;
    }
  }
}


template<typename valtype, typename indtype, bool mk> // this function is for experimenting
inline valtype **mvalUpperBoundBiMan(valtype **begin, valtype **end, valtype *val,
                                     indtype compst, indtype dcomp, INT *mask)
// find the first one that is notAllLessEqual
{
  --end;
  std::swap(begin, end);
  // if(!notAllLessEqual(*begin, val, d)) return begin + 1;
  if(!notAllLessEqual<valtype, indtype, mk> (
      *begin + compst, val + compst, dcomp, mask)) return begin + 1;
  while(true)
  {
    valtype **mid = begin - int((begin - end) / 2);
    // if(notAllLessEqual(*mid, val, d))
    if(notAllLessEqual<valtype, indtype, mk> (
        *mid + compst, val + compst, dcomp, mask))
    {
      if(mid == begin)
      {
        return end + 1;
      }
      begin = mid;
    }
    // to here, it means mid is not all less equal than val
    // else if(!notAllLessEqual(*(mid + 1), val, d))
    else if(!notAllLessEqual<valtype, indtype, mk> (
        *(mid + 1) + compst, val + compst, dcomp, mask))
    // and at the same time even mid-1 is not all less equal than val
    {
      end = mid;// so end is safe to set as mid
    }
    // to here, it means mid is the first that not all less equal than val from low to high
    else
    {
      return mid + 1;
    }
  }
}




template<typename valtype, typename indtype, bool mk>
inline valtype **mvalLowerBoundLr(valtype **begin, valtype **end, valtype *val,
                                  indtype compst, indtype dcomp, INT *mask)
{
  indtype ic = 0;
  for(; begin < end; ++begin)
  {
    if(!notAllGreaterEqual<valtype, indtype, mk> (
        *begin + compst, val + compst, ic, dcomp, mask))
    {
      break;
    }
  }
  return begin;
}


template<typename valtype, typename indtype, bool mk>
inline valtype **mvalLowerBoundLr(valtype **begin, valtype **end, valtype *val, indtype &ic,
                                  indtype compst, indtype dcomp, INT *mask) // ic is relative to compst
{
  for(; begin < end; ++begin)
  {
    if(!notAllGreaterEqual<valtype, indtype, mk> (
        *begin + compst, val + compst, ic, dcomp, mask))
    {
      break;
    }
  }
  return begin;
}


template<typename valtype, typename indtype, bool mk>
inline valtype **mvalUpperBoundLr(valtype **begin, valtype **end, valtype *val,
                                  indtype compst, indtype dcomp, INT *mask)
{
  indtype ic = 0;
  valtype **i = end - 1;
  for(; i >= begin; --i)
  {
    if(!notAllLessEqual<valtype, indtype, mk> (
        *i + compst, val + compst, ic, dcomp, mask))
    {
      break;
    }
  }
  return i;
}


template<typename valtype, typename indtype, bool mk>
inline valtype **mvalUpperBoundLr(valtype **begin, valtype **end, valtype *val, indtype &ic,
                                  indtype compst, indtype dcomp, INT *mask)
{
  valtype **i = end - 1;
  for(; i >= begin; --i)
  {
    if(!notAllLessEqual<valtype, indtype, mk> (
        *i + compst, val + compst, ic, dcomp, mask))
    {
      break;
    }
  }
  return i;
}



