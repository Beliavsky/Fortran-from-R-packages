#pragma once
#include <Rcpp.h>
// using namespace Rcpp;


template<typename ing>
inline void raiseSeqBy1(ing *begin, ing *end)
{
  for(; begin < end; ++begin)
  {
    ++(*begin);
  }
}


template<typename ing>
inline void supressSeqBy1(ing *begin, ing *end)
{
  for(; begin < end; ++begin)
  {
    --(*begin);
  }
}


template<typename ing>
inline ing riseUpdateRightConsecutiveSeqEnd(
    ing *base, ing *cnsctvSeqEnd, ing *seqRightEnd)
// Find the consecutive sequence's end on the right side of LB,
// specifically for position = 0 in the update procedure
// and the right part of LB when position is in middle
{
  while(cnsctvSeqEnd < seqRightEnd)
  {
    if(*cnsctvSeqEnd == cnsctvSeqEnd[-1] + 1)
    {
      ++cnsctvSeqEnd;
    }
    else break;
  }
  return cnsctvSeqEnd - base;
}


template<typename ing>
inline ing descendUpdateLeftConsecutiveSeqEnd(
    ing *base, ing *cnsctvSeqEnd, ing *seqLeftEnd)
// Find the consecutive sequence's end on the left side of UB, specifically for
// position = len in the update procedure.
{
  while(cnsctvSeqEnd > seqLeftEnd)
  {
    if(*cnsctvSeqEnd == cnsctvSeqEnd[1] - 1)
    {
      --cnsctvSeqEnd;
    }
    else break;
  }
  return cnsctvSeqEnd - base;
}


template<typename ing>
inline ing riseUpdateLeftConsecutiveSeqBegin(
    ing *base, ing *cnsctvSeqBegin,
    ing *seqRightEnd, ing *UBleftResvBegin)
  // Find the consecutive sequence's begin on the left side of UB, specifically
  // for position is in middle in the update()
{
  ing *tmp = UBleftResvBegin + (cnsctvSeqBegin - base);
  while(cnsctvSeqBegin < seqRightEnd)
  {
    if(*cnsctvSeqBegin == *tmp)
    {
      ++tmp;
      ++cnsctvSeqBegin;
    }
    else break;
  }
  return cnsctvSeqBegin - base;
}


template<typename ing>
inline ing consecutiveBoundSeqFromRightToLeft(
    ing *base, ing *seqRightBegin, ing boundRightBegin)
{
  if(*seqRightBegin <= boundRightBegin)
  {
    return seqRightBegin + 1 - base;
  }


  *seqRightBegin = boundRightBegin;
  while(seqRightBegin > base)
  {
    // if(seqRightBegin <= base) break;
    --boundRightBegin;
    if(seqRightBegin[-1] > boundRightBegin)
    {
      --seqRightBegin;
      *seqRightBegin = boundRightBegin;
    }
    else break;
  }
  return seqRightBegin - base;
}










