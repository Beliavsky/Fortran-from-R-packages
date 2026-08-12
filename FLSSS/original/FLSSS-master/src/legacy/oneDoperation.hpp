# pragma once
// # include <Rcpp.h>
# include "macros.hpp"
// using namespace Rcpp;


template<typename valtype, typename indtype>
inline void iterSum(valtype &rst, valtype *V, indtype *ele, indtype eleSize)
{
  for(indtype i = 0; i < eleSize; ++i)
  {
    rst += V[ele[i]];
  }
}


template<typename valtype>
inline valtype absDiff(valtype x, valtype y)
{
  if(y == 0)
  {
    if(x == 0) return 0;
    return infy;
  }
  return std::abs(x / y - 1);
}


//—————————————————————————————————————————————————————————————————————————————————————————————————
template<typename valtype, typename indtype> // this function is for experimenting
inline valtype *lowerBoundBiMan(valtype *begin, valtype *end, valtype val)
{
  if(*begin >= val) return begin;
  --end;
  while(true)
  {
    valtype *mid = begin + int((end - begin) / 2);
    // if(notAllGreaterEqual(*mid, val, d))
    if(*mid < val)
    {
      if(mid == begin)
      {
        return end;
      }
      begin = mid;
    }
    // to here, it means mid is all greater equal than val
    else if(*(mid - 1) >= val)
    // and at the same time even mid-1 is all greater equal than val
    {
      end = mid;// so end is safe to set as mid
    }
    // to here, it means mid is the first that all greater equal than val from low to high
    else
    {
      return mid;
    }
  }
}


template<typename valtype, typename indtype> // this function is for experimenting
inline valtype *upperBoundBiMan(valtype *begin, valtype *end, valtype val)
// find the first one that is notAllLessEqual
{
  --end;
  std::swap(begin, end);
  if(*begin <= val) return begin + 1;
  while(true)
  {
    valtype *mid = begin - int((begin - end) / 2);
    if(*mid > val)
    {
      if(mid == begin)
      {
        return end + 1;
      }
      begin = mid;
    }
    // to here, it means mid is not all less equal than val
    else if(*(mid + 1) <= val)
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
inline valtype *lowerBoundLr(valtype *begin, valtype *end, valtype val)
{
  for(; begin < end; ++begin)
  {
    if(*begin >= val)
    {
      break;
    }
  }
  return begin;
}


template<typename valtype, typename indtype>
inline valtype *upperBoundLr(valtype *begin, valtype *end, valtype val)
{
  for(--end; end >= begin; --end)
  {
    if(*end <= val)
    {
      break;
    }
  }
  return end;
}



