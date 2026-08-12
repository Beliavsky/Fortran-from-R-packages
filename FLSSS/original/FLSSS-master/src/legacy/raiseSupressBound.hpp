# pragma once
# include <Rcpp.h>


namespace legacy
{
template<typename indtype>
inline void raiseSeqBy1(indtype *begin, indtype *end)
{
  for(; begin < end; ++begin)
  {
    ++(*begin);
  }
}


template<typename indtype>
inline void supressSeqBy1(indtype *begin, indtype *end)
{
  for(; begin < end; ++begin)
  {
    --(*begin);
  }
}


template<typename indtype>
inline indtype riseUpdateRightConsecutiveSeqEnd(indtype *base, indtype *cnsctvSeqEnd, indtype *seqRightEnd)
// find the consecutive sequence's end on the right side of LB, specifically for position = 0 in the update procedure
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


template<typename indtype>
inline indtype descendUpdateLeftConsecutiveSeqEnd(indtype *base, indtype *cnsctvSeqEnd, indtype *seqLeftEnd)
// find the consecutive sequence's end on the left side of UB, specifically for position=len in the update procedure
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


template<typename indtype>
inline indtype riseUpdateLeftConsecutiveSeqBegin(indtype *base, indtype *cnsctvSeqBegin,
                                                 indtype *seqRightEnd, indtype *UBleftResvBegin)
  // find the consecutive sequence's begin on the left side of UB, specifically for position is in middle in the update()
{
  indtype *tmp = UBleftResvBegin + (cnsctvSeqBegin - base);
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


template<typename indtype>
inline indtype consecutiveBoundSeqFromRightToLeft(indtype *base, indtype *seqRightBegin, indtype boundRightBegin)
{
  if(*seqRightBegin <= boundRightBegin)
  {
    return seqRightBegin + 1 - base;
  }


  *seqRightBegin = boundRightBegin;
  while(true)
  {
    if(seqRightBegin <= base) break;
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


}

