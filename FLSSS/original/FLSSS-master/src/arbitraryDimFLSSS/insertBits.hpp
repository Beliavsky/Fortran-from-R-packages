#pragma once


// Insert at most 64 bits in y into x.
// Return true if y has been exhausted.
inline void insertBits(uint64_t *x, int &XbeginBit,
                       uint64_t *y, int &YbeginBit,
                       int YbitEndWhichInt, int YbitEndWhichBit)
{
  constexpr int NbitPerInt = sizeof(uint64_t) * 8;
  int whichIntX = XbeginBit / NbitPerInt,
    whichBitX = XbeginBit % NbitPerInt,
    whichIntY = YbeginBit / NbitPerInt,
    whichBitY = YbeginBit % NbitPerInt;
  uint64_t addition = y[whichIntY];
  int NbitsFilled = 0;
  if(YbitEndWhichInt == whichIntY) // Bits in Y end in this integer.
  {
    int tmp = NbitPerInt - YbitEndWhichBit;
    addition = (addition << tmp) >> tmp;
    addition = (addition >> whichBitY) << whichBitX;
    x[whichIntX] += addition;
    NbitsFilled = std::min(NbitPerInt - whichBitX, YbitEndWhichBit - whichBitY);
  }
  else
  {
    addition = (addition >> whichBitY) << whichBitX;
    x[whichIntX] += addition;
    NbitsFilled = std::min(NbitPerInt - whichBitX, NbitPerInt - whichBitY);
  }
  XbeginBit += NbitsFilled;
  YbeginBit += NbitsFilled;
}


// NeffectiveBits[i] is the number of effective bits in the array of y[i].
// x[] should have been zeroed, and will contain the result.
inline void insertAllBits(uint64_t *x, vec<vec<uint64_t> > &y,
                          int *NeffectiveBits)
{
  int XbeginBit = 0;
  for(int i = 0, iend = y.size(); i < iend; ++i)
  {
    int YbeginBit = 0, XbeginBitLast = XbeginBit;
    int YendBitWhichInt = NeffectiveBits[i] / (sizeof(uint64_t) * 8);
    int YendBitWhichBit = NeffectiveBits[i] % (sizeof(uint64_t) * 8);
    while(XbeginBit - XbeginBitLast < NeffectiveBits[i])
      insertBits(x, XbeginBit, &y[i][0], YbeginBit,
                 YendBitWhichInt, YendBitWhichBit);
  }
}
// Test kit in ../arbitraryDimFLSSStest/insertBits.cpp








