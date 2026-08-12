#pragma once
#include "charlieThreadPool.hpp"
#ifndef vec
#define vec std::vector
#endif
#include "arithmetic.hpp"
#include "insertBits.hpp"
#include "validateStringInput.hpp"


template<typename T>
inline void removeHighEndZeros(vec<T> &x)
{
  int i = x.size() - 1;
  for(; i >= 1 and x[i] == 0; --i);
  x.resize(i + 1); // x should be of size 1 at least, if x is nonempty.
}


// x's size can be adjusted.
template<typename T>
inline void normalizeDecimalIntVec(vec<T> &v)
{
  for(int k = 0, kend = v.size() - 1; k < kend; ++k)
  {
    if(v[k] >= 0)
    {
      v[k + 1] += v[k] / 10;
      v[k] %= 10;
    }
    else
    {
      T subtractor = (-v[k] + 9) / 10;
      v[k + 1] -= subtractor;
      v[k] += subtractor * 10;
    }
  }
  while(v.back() > 9)
  {
    int vk = v.back();
    v.push_back(vk / 10);
    v[v.size() - 2] = vk % 10;
  }
  removeHighEndZeros<T> (v);
}


template<typename T>
inline bool nonnegVecLessThan(T *x, int xsize, T *y, int ysize)
{
  if(xsize < ysize) return true;
  if(xsize > ysize) return false;
  for(int k = xsize - 1; k >= 0; --k)
  {
    if(x[k] < y[k]) return true;
    if(x[k] > y[k]) return false;
  }
  return false;
}


template<typename T>
struct ComparePosiVec
{
  vec<T> *val;
  ComparePosiVec(vec<T> *val): val(val) {}
  bool operator()(int i, int j)
  {
    return nonnegVecLessThan<T> (
      &val[i][0], val[i].size(), &val[j][0], val[j].size());
  }
};


struct AdjustStringVec
{
  vec<char> sign;
  vec<int> Ndeci;
  vec<vec<int> > value;
  vec<int> targetSum;


  // return x < y
  bool valuevecLessThan(int *x, int xsize, int *y,
                        int ysize, char xsign, char ysign)
  {
    if(xsign == -1 and ysign == 1) return true;
    if(xsign == 1 and ysign == -1) return false;
    if(xsign == -1 and ysign == -1)
    {
      std::swap(x, y);
      std::swap(xsize, ysize);
    }
    if(xsize < ysize) return true;
    if(xsize > ysize) return false;
    int i = xsize - 1;
    for(; i >= 0; --i)
    {
      if(x[i] < y[i]) return true;
      if(x[i] > y[i]) return false;
    }
    return false;
  }


  // len is the subset size.
  vec<char*> x; // Pointers to value strings. Last element should be the target sum.
  vec<int> xsizes; // Sizes of those strings.
  // Xval and tsum are references.
  // vec<unsigned char> valid; // true if the element is <= target sum.


  void stringsToBitIntVecInitialize(StringMatrix &Xval, int whichCol,
                                    StringVector &tsum, int len)
  {
    x.resize(Xval.nrow() + 1); // The last one is to store tsum.
    xsizes.resize(x.size());
    StringMatrix::Column xval = Xval(_, whichCol);
    for (int i = 0, iend = x.size() - 1; i < iend; ++i)
    {
      const char *a = &xval[i][0];
      std::memcpy(&x[i], &a, sizeof(char*));
      xsizes[i] = xval[i].size();
    }
    const char *a = &tsum[whichCol][0];
    std::memcpy(&x.back(), &a, sizeof(char*));
    xsizes.back() = tsum[whichCol].size();
  }


  void stringsToBitIntVec(int len)
  {
    sign.assign(x.size(), 1);
    Ndeci.resize(x.size());
    int maxDecimal = 0;
    for (int i = 0, iend = x.size(); i < iend; ++i)
    {
      char *u = x[i];
      int usize = xsizes[i];
      sign[i] = (u[0] != '-') * 2 - 1;
      int k = 0;
      for (; k < usize and u[k] != '.'; ++k);
      Ndeci[i] = std::max(usize - 1 - k, 0);
      maxDecimal = std::max<int> (Ndeci[i], maxDecimal);
    }


    vec<int> &zerosToAdd = Ndeci;
    for (int i = 0, iend = zerosToAdd.size(); i < iend; ++i)
    {
      zerosToAdd[i] = maxDecimal - zerosToAdd[i];
    }


    value.resize(zerosToAdd.size());
    for (int i = 0, iend = zerosToAdd.size(); i < iend; ++i)
    {
      char *u = x[i];
      int usize = xsizes[i];
      value[i].reserve(usize + zerosToAdd[i]);
      value[i].assign(zerosToAdd[i], 0);
      for (int k = usize - 1; k >= 0; --k)
      {
        if (u[k] != '-' and u[k] != '.') value[i].push_back(u[k] - 48);
      }
    }


    // Find minimum.
    int mini = 0;
    for (int i = 1, iend = sign.size() - 1; i < iend; ++i) // Do not include target.
    {
      vec<int> &minx = value[mini];
      vec<int> &current = value[i];
      bool currentIsLess = valuevecLessThan(
        &current[0], current.size(), &minx[0], minx.size(), sign[i], sign[mini]);
      if (currentIsLess) mini = i;
    }


    // Shift value.
    if(sign[mini] == -1)
    {
      int *u = &value[mini][0];
      for(int k = 0, kend = value[mini].size(); k < kend; ++k) u[k] = -u[k];
    }
    for(int i = 0, iend = sign.size(); i < iend; ++i)
    {
      if(i == mini) continue;
      int multiple = i == iend - 1 ? len : 1; // For adjusting the target sum.
      vec<int> &v = value[i];
      vec<int> &vmin = value[mini];
      v.resize(std::max(vmin.size(), v.size()));
      if(sign[i] == -1)
      {
        for(int k = 0, kend = v.size(); k < kend; ++k)
        {
          v[k] = -v[k];
        }
      }
      for(int k = 0, kend = std::min(v.size(), vmin.size()); k < kend; ++k)
      {
        v[k] -= vmin[k] * multiple;
      }
      normalizeDecimalIntVec<int> (v);
    }
    value[mini].assign(1, 0);


    targetSum.swap(value.back());
    value.resize(value.size() - 1);
  }


  void subtractMinAdjustSum(int len) // Only after all entries are nonnegative.
  {
    // Find minimum.
    int mini = 0;
    for(int i = 1, iend = value.size(); i < iend; ++i)
    {
      vec<int> &minx = value[mini];
      vec<int> &current = value[i];
      bool currentIsLess = valuevecLessThan(
        &current[0], current.size(), &minx[0], minx.size(), 1, 1);
      if(currentIsLess) mini = i;
    }


    for(int i = 0, iend = value.size(); i < iend; ++i)
    {
      if(i == mini) continue;
      vec<int> &v = value[i], &vmin = value[mini];
      v.resize(std::max(vmin.size(), v.size()));
      for(int k = 0, kend = std::min(v.size(), vmin.size()); k < kend; ++k)
        v[k] -= vmin[k];
      normalizeDecimalIntVec<int> (v);
    }


    targetSum.resize(std::max(value[mini].size(), targetSum.size()));
    for(int k = 0, kend = std::min(targetSum.size(), value[mini].size()); k < kend; ++k)
      targetSum[k] -= value[mini][k] * len;
    normalizeDecimalIntVec<int> (targetSum);


    value[mini].assign(1, 0);
  }


  vec<int> rank, order;
  void createRankOrder()
  {
    order.resize(value.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), ComparePosiVec<int> (&value[0]));
    rank.resize(order.size());
    for(int i = 0, iend = rank.size(); i < iend; ++i)
      rank[order[i]] = i;
  }


  vec<int> largestSubsetSum;


  // USE THIS FUNCTION IF K-SUM IS MADE INTERNALLY.
  void getLargestSubsetSum(int subsetlen)
  {
    int veclen = 0;
    for(int i = order.size() - 1, iend = order.size() - subsetlen; i >= iend; --i)
      veclen = std::max<int> (veclen, value[order[i]].size());
    largestSubsetSum.resize(veclen, 0);


    for(int i = order.size() - 1, iend = order.size() - subsetlen; i >= iend; --i)
    {
      int k = order[i];
      for(int j = 0, jend = value[k].size(); j < jend; ++j)
      {
        largestSubsetSum[j] += value[k][j];
      }
    }
    normalizeDecimalIntVec(largestSubsetSum);
  }


  // USE THIS FUNCTION IS K-SUM IS IMPORTED OR TO BE EXPORTED.
  void getLargestSubsetSum()
  {
    int veclen = 0;
    for (int i = 0, iend = order.size(); i < iend; ++i)
      veclen = std::max<int> (veclen, value[i].size());
    largestSubsetSum.resize(veclen, 0);
    for(int i = 0, iend = order.size(); i < iend; ++i)
    {
      for(int j = 0, jend = value[i].size(); j < jend; ++j)
        largestSubsetSum[j] += value[i][j];
    }
    normalizeDecimalIntVec(largestSubsetSum);
  }


  void bitIntsToHugeInts64(vec<uint64_t> &rst, vec<int> &binIntVec)
  {
    int ssize = binIntVec.size();
    if(ssize <= 0) return;
    int *s = &binIntVec[0];
    int Nbit = int(ssize * 3.33) + 1;
    int rstSize = (Nbit + 31) / 32; // The number of 32-bit integers.
    rst.resize((rstSize + 1) / 2, 0);
    vec<uint64_t> buffer((rstSize * 3 + 1) / 2, 0);
    uint32_t *power10 = (uint32_t*)(&buffer[0]),
      *power10new = power10 + rstSize, *adder = power10new + rstSize;
    power10[0] = 1;
    uint32_t ten = 10;
    // for(int i = ssize - 1; ; --i)
    for(int i = 0; ; ++i)
    {
      uint32_t a = s[i];
      mulHint32(adder, rstSize, &a, 1, power10, rstSize);
      addHint32((uint32_t*)&rst[0], (uint32_t*)&rst[0], adder, rstSize);
      // if(i <= 0) break;
      if(i >= ssize - 1) break;
      mulHint32(power10new, rstSize, power10, rstSize, &ten, 1);
      std::swap(power10new, power10);
    }
    removeHighEndZeros<uint64_t> (rst);
  }


  vec<vec<uint64_t> > value64bit;
  vec<uint64_t> targetSum64bit;
  vec<uint64_t> largestSubsetSum64bit;
  int generate64bitUint() // Return the number of bits in largest subset sum.
  {
    bitIntsToHugeInts64(largestSubsetSum64bit, largestSubsetSum);


    value64bit.resize(value.size());
    for(int i = 0, iend = value.size(); i < iend; ++i)
    {
      bitIntsToHugeInts64(value64bit[i], value[i]);
      while(value64bit[i].size() < largestSubsetSum64bit.size())
        value64bit[i].push_back(0);
    }


    bitIntsToHugeInts64(targetSum64bit, targetSum);
    while(targetSum64bit.size() < largestSubsetSum64bit.size())
      targetSum64bit.push_back(0);


    int uselessHighBits = 0;
    if(true)
    {
      uint64_t &x = largestSubsetSum64bit.back();
      for(int shift = sizeof(uint64_t) * 8 - 1; shift >= 0; --shift)
      {
        if((x >> shift) != 0) break;
        uselessHighBits += 1;
      }
    }
    int validbits = sizeof(uint64_t) * 8 * largestSubsetSum64bit.size() - uselessHighBits;
    return validbits;
  }


};




// Return true if the matrix is not empty.
// order is the row order.
inline bool stringMatTo64bitIntMat(
    StringMatrix &Xval, StringVector &tsum, int len, CharlieThreadPool &tp,
    vec<uint64_t> &X, vec<uint64_t> &targetSS, vec<uint64_t> &largestSS,
    vec<int> &order, vec<int> &colOrder, int &totalEffectiveBits,
    bool ksumIsComputedInternallyAndWillNotBeExported)
{


  validateNumStrings(Xval);
  validateNumStrings(tsum);


  auto empf = [](std::size_t t){ return false; };


  vec<AdjustStringVec> adjv(tsum.size());
  for(int i = 0, iend = adjv.size(); i < iend; ++i)
    adjv[i].stringsToBitIntVecInitialize(Xval, i, tsum, len);


  tp.parFor(0, adjv.size(), [&](std::size_t i, std::size_t t)
  {
    adjv[i].stringsToBitIntVec(len);
    return false;
  }, 1, empf, empf);


  // // If targetSum of a certain dimension is negative, solution does not exist.
  bool rstbool = true;
  for(int i = 0, iend = adjv.size(); i < iend; ++i)
  {
    if(adjv[i].targetSum.back() < 0)
    {
      rstbool = false;
      break; // DO NOT RETURN IMMEDIATELY. THE REST IS NEEDED FOR K-SUM!
    }
  }
  // Yes, the above is true, but for ksum, we need the bit matrix of superset
  // any way, so stick to it.


  // Compute rank orders and largest subset sums.
  tp.parFor(0, adjv.size(), [&](std::size_t i, std::size_t t)
  {
    adjv[i].createRankOrder();
    if (ksumIsComputedInternallyAndWillNotBeExported)
      adjv[i].getLargestSubsetSum(len); // adjv[i].getLargestSubsetSum();
    else adjv[i].getLargestSubsetSum();
    return false;
  }, 1, empf, empf);


  // Compute rank correlations with others (just inner product)
  vec<std::size_t> innerProd(adjv.size(), 0);
  tp.parFor(0, adjv.size(), [&](std::size_t i, std::size_t t)
  {
    auto &x = adjv[i];
    for(int k = 0, kend = adjv.size(); k < kend; ++k)
      innerProd[i] += std::inner_product(
        x.rank.begin(),  x.rank.end(), adjv[k].rank.begin(), 0);
    return false;
  }, 1, empf, empf);


  // vec<int> colOrder(adjv.size());
  colOrder.resize(adjv.size());
  std::iota(colOrder.begin(), colOrder.end(), 0);
  std::sort(colOrder.begin(), colOrder.end(),
            [&](const int &x, const int &y)->bool
  {
    return innerProd[x] < innerProd[y];
  });


  // Generate 64bit integers and then reorder the value vectors.
  vec<int> Nbits(adjv.size());
  tp.parFor(0, adjv.size(), [&](std::size_t i, std::size_t t)
  {
    Nbits[i] = adjv[i].generate64bitUint();
    return false;
  }, 1, empf, empf);


  // # of 64-bit integers in every row.
  int totalBits = std::accumulate(Nbits.begin(), Nbits.end(), 0);
  uint64_t dim = (totalBits + sizeof(uint64_t) * 8 - 1) / (sizeof(uint64_t) * 8);
  if(totalBits % (sizeof(uint64_t) * 8) == 0)
    // Unfortunately there needs to be one more integer for storing the sign bit.
  {
    totalBits += 1;
    dim += 1;
  }


  totalEffectiveBits = totalBits;
  X.resize(adjv[0].value64bit.size() * dim);


  vec<int> NbitsColOrdered(Nbits.size());
  for(int i = 0, iend = Nbits.size(); i < iend; ++i)
    NbitsColOrdered[i] = Nbits[colOrder[i]];
  Nbits.swap(NbitsColOrdered);


  // Put all the bits in X.
  tp.parFor(0, adjv[0].value.size(), [&](std::size_t i, std::size_t t)
  {
    vec<vec<uint64_t> > row(adjv.size());
    for(int j = 0, jend = row.size(); j < jend; ++j)
      row[j].swap(adjv[ colOrder[j] ].value64bit[i]);
    insertAllBits(&X[0] + dim * i, row, &Nbits[0]);
    return false;
  }, 1, empf, empf);


  targetSS.resize(dim);
  largestSS.resize(dim);
  if(true)
  {
    vec<vec<uint64_t> > tmp(adjv.size());


    for(int i = 0, iend = adjv.size(); i < iend; ++i)
      tmp[i].swap(adjv[ colOrder[i] ].targetSum64bit);
    insertAllBits(&targetSS[0], tmp, &Nbits[0]);


    for(int i = 0, iend = adjv.size(); i < iend; ++i)
      tmp[i].swap(adjv[ colOrder[i] ].largestSubsetSum64bit);
    insertAllBits(&largestSS[0], tmp, &Nbits[0]);
  }


  order.resize(X.size() / dim);
  std::iota(order.begin(), order.end(), 0);
  sort(order.begin(), order.end(), [&](const int &x, const int &y)->bool
  {
    uint64_t *a = &X[0] + x * dim, *b = &X[0] + y * dim;
    for(int i = dim - 1; i >= 0; --i)
    {
      if(a[i] < b[i]) return true;
      if(a[i] > b[i]) return false;
    }
    return false;
  });


  vec<uint64_t> Xtmp(X.size());
  for(int i = 0, iend = order.size(); i < iend; ++i)
    std::memcpy(&Xtmp[0] + dim * i, &X[0] + dim * order[i],
                sizeof(uint64_t) * dim);
  X.swap(Xtmp);


  return rstbool;
}




































































































