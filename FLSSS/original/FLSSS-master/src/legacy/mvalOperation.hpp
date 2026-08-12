# pragma once
// # include <Rcpp.h>
# include "macros.hpp"
// using namespace Rcpp;


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
    *i = *left - *right + *mid;
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


template<typename valtype, typename indtype>
inline bool notAllLessEqual(valtype *x, valtype *y, indtype &i, indtype d)
{
  for(; i < d; ++i)
  {
    if(x[i] > y[i]) return 1;
  }
  return 0;
}


template<typename valtype, typename indtype>
inline bool notAllLessEqual(valtype *x, valtype *y, indtype d)
{
  for(indtype i = 0; i < d; ++i)
  {
    if(x[i] > y[i]) return 1;
  }
  return 0;
}


template<typename valtype, typename indtype>
inline bool notAllGreaterEqual(valtype *x, valtype *y, indtype &i, indtype d)
{
  for(; i < d; ++i)
  {
    if(x[i] < y[i]) return 1;
  }
  return 0;
}


template<typename valtype, typename indtype>
inline bool notAllGreaterEqual(valtype *x, valtype *y, indtype d)
{
  for(indtype i = 0; i < d; ++i)
  {
    if(x[i] < y[i]) return 1;
  }
  return 0;
}


template<typename valtype, typename indtype>
inline bool allEqual(valtype *x, valtype *y, indtype &i, indtype d)
{
  for(; i < d; ++i)
  {
    if(std::abs(x[i] - y[i]) > eps) return 0;
  }
  return 1;
}


template<typename valtype, typename indtype>
inline bool allEqual(valtype *x, valtype *y, indtype d)
{
  for(indtype i = 0; i < d; ++i)
  {
    if(std::abs(x[i] - y[i]) > eps) return 0;
  }
  return 1;
}








// ————————————————————————————————————————————————————————————————————————————————————————————————
/*
template<typename valtype, typename indtype>
struct mvalLowerThan
{
  indtype &d;
  mvalLowerThan(indtype &d):d(d){};
  bool operator() (valtype *u, valtype *v)
  {
    for(valtype *uend = u + d; u < uend; ++u, ++v)
    {
      if(*u < *v) return 1; // as long as there's one less, then jump up
    }
    return 0;
  }
};


template<typename valtype, typename indtype>
inline valtype **mvalLowerBoundBi(valtype **begin, valtype **end, valtype *val, indtype d)
{
  return std::lower_bound(begin, end, val, mvalLowerThan <valtype,indtype> (d));
}


template<typename valtype, typename indtype>
inline valtype **mvalUpperBoundBi(valtype **begin, valtype **end, valtype *val, indtype d)
{
  return std::upper_bound(begin, end, val, mvalLowerThan<valtype, indtype> (d));
}
*/




template<typename valtype, typename indtype> // this function is for experimenting
// inline valtype **mvalLowerBoundBiMan(valtype **begin, valtype **end, valtype *val, indtype d)
inline valtype **mvalLowerBoundBiMan(valtype **begin, valtype **end, valtype *val,
                                     indtype compst, indtype compsiz)
{
  // if(!notAllGreaterEqual(*begin, val, d)) return begin;
  if(!notAllGreaterEqual(*begin + compst, val + compst, compsiz)) return begin;
  --end;
  while(true)
  {
    valtype **mid = begin + int((end - begin) / 2);
    // if(notAllGreaterEqual(*mid, val, d))
    if(notAllGreaterEqual(*mid + compst, val + compst, compsiz))
    {
      if(mid == begin)
      {
        return end;
      }
      begin = mid;
    }
    // to here, it means mid is all greater equal than val
    // else if(!notAllGreaterEqual(*(mid - 1), val, d))
    else if(!notAllGreaterEqual(*(mid - 1) + compst, val + compst, compst))
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


template<typename valtype, typename indtype> // this function is for experimenting
inline valtype **mvalUpperBoundBiMan(valtype **begin, valtype **end, valtype *val,
                                     indtype compst, indtype dcomp)
// find the first one that is notAllLessEqual
{
  --end;
  std::swap(begin, end);
  // if(!notAllLessEqual(*begin, val, d)) return begin + 1;
  if(!notAllLessEqual(*begin + compst, val + compst, dcomp)) return begin + 1;
  while(true)
  {
    valtype **mid = begin - int((begin - end) / 2);
    // if(notAllLessEqual(*mid, val, d))
    if(notAllLessEqual(*mid + compst, val + compst, dcomp))
    {
      if(mid == begin)
      {
        return end + 1;
      }
      begin = mid;
    }
    // to here, it means mid is not all less equal than val
    // else if(!notAllLessEqual(*(mid + 1), val, d))
    else if(!notAllLessEqual(*(mid + 1) + compst, val + compst, dcomp))
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




template<typename valtype, typename indtype>
inline valtype **mvalLowerBoundLr(valtype **begin, valtype **end, valtype *val,
                                  indtype compst, indtype dcomp)
{
  indtype ic = 0;
  for(; begin < end; ++begin)
  {
    // if(!notAllGreaterEqual(*begin, val, ic, d))
    if(!notAllGreaterEqual(*begin + compst, val + compst, ic, dcomp))
    {
      break;
    }
  }
  return begin;
}


template<typename valtype, typename indtype>
inline valtype **mvalLowerBoundLr(valtype **begin, valtype **end, valtype *val, indtype &ic,
                                  indtype compst, indtype dcomp) // ic is relative to compst
{
  for(; begin < end; ++begin)
  {
    if(!notAllGreaterEqual(*begin + compst, val + compst, ic, dcomp))
    {
      break;
    }
  }
  return begin;
}


template<typename valtype, typename indtype>
inline valtype **mvalUpperBoundLr(valtype **begin, valtype **end, valtype *val,
                                  indtype compst, indtype dcomp)
{
  indtype ic = 0;
  valtype **i = end - 1;
  for(; i >= begin; --i)
  {
    if(!notAllLessEqual(*i + compst, val + compst, ic, dcomp))
    {
      break;
    }
  }
  return i;
}


template<typename valtype, typename indtype>
inline valtype **mvalUpperBoundLr(valtype **begin, valtype **end, valtype *val, indtype &ic,
                                  indtype compst, indtype dcomp)
{
  valtype **i = end - 1;
  for(; i >= begin; --i)
  {
    if(!notAllLessEqual(*i + compst, val + compst, ic, dcomp))
    {
      break;
    }
  }
  return i;
}



