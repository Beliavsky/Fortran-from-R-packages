// [[Rcpp::plugins(cpp17)]]
# include <Rcpp.h>
# include <bitset>
// [[Rcpp::depends(RcppParallel)]]
# include "../arbitraryDimFLSSS/dnyTasking2.hpp"
# include "../arbitraryDimFLSSS/arithmetic.hpp"
using namespace Rcpp;
# include "../arbitraryDimFLSSS/validateStringInput.hpp"




#ifndef vec
#define vec std::vector
#endif


// z can equal x or y.
void addHint32(uint32_t *z, uint32_t *x, uint32_t *y, int size)
{
  uint64_t s = uint64_t(x[0]) + y[0];
  z[0] = uint32_t(s);
  for(int i = 1; i < size; ++i)
  {
    s = (s >> 32) + x[i] + y[i];
    z[i] = uint32_t(s);
  }
}


// z cannot equal x or y.
// z does not have be initialized by zeros.
void mulHint32(uint32_t *z, int zsize,
               uint32_t *x, int xsize,
               uint32_t *y, int ysize)
{
  for(--xsize; xsize >= 0 and x[xsize] == 0; --xsize);
  ++xsize;
  for(--ysize; ysize >= 0 and y[ysize] == 0; --ysize);
  ++ysize;
  std::fill(z, z + zsize, 0);
  for(int i = 0; i < xsize; ++i)
  {
    for(int j = 0; j < ysize; ++j)
    {
      int k = i + j;
      uint64_t s = uint64_t(x[i]) * y[j] + z[k];
      z[k] = uint32_t(s);
      s >>= 32;
      for(++k; s != 0; ++k, s >>= 32) // s is the carry;
      {
        s += z[k];
        z[k] = uint32_t(s);
      }
    }
  }
}


// [[Rcpp::export]]
List testMulHintVbinInput(IntegerVector xv, IntegerVector yv)
{
  vec<int> x(xv.begin(), xv.end()), y(yv.begin(), yv.end());
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());
  int asize = (x.size() + 31) / 32;
  int bsize = (y.size() + 31) / 32;
  vec<uint32_t> a(asize, 0), b(bsize, 0);
  for(int i = 0, k = 0; i < asize; ++i)
  {
    for(int j = 0; j < 32; ++k, ++j)
      a[i] += ((uint32_t)x[k]) << j;
  }
  for(int i = 0, k = 0; i < bsize; ++i)
  {
    for(int j = 0; j < 32; ++k, ++j)
      b[i] += ((uint32_t)y[k]) << j;
  }


  if(false)
  {
    std::string astr = "";
    for(int i = a.size() - 1; i >= 0; --i)
      astr += std::bitset<32>(a[i]).to_string();
    std::string bstr = "";
    for(int i = b.size() - 1; i >= 0; --i)
      bstr += std::bitset<32>(b[i]).to_string();
    std::cout << "astr = " << astr << std::endl;
    std::cout << "bstr = " << bstr << std::endl;
  }


  vec<uint32_t> ab(a.size() + b.size(), 0);
  // return List::create();
  mulHint32(&ab[0], ab.size(), &a[0], a.size(), &b[0], b.size());
  std::string abstr = "";
  for(int i = ab.size() - 1; i >= 0; --i)
    abstr += std::bitset<32>(ab[i]).to_string();


  vec<int> conv(x.size() + y.size(), 0);
  for(int i = 0, iend = x.size(); i < iend; ++i)
  {
    for(int j = 0, jend = y.size(); j < jend; ++j)
    {
      conv[i + j] += x[i] * y[j];
    }
  }
  for(int i = 0, iend = conv.size() - 1; i < iend; ++i)
  {
    conv[i + 1] += conv[i] / 2;
    conv[i] %= 2;
  }


  std::reverse(conv.begin(), conv.end());


  std::string correctRst(conv.size(), '0');
  for(int i = 0, iend = conv.size(); i < iend; ++i)
    correctRst[i] = conv[i] == 1 ? '1' : '0';


  return List::create(Named("rst") = abstr, Named("correctRst") = correctRst);
}


void intString2hugeInt32(vec<uint32_t> &rst, char *s, int ssize)
{
  if(ssize <= 0) return;
  int Nbit = int(ssize * 3.33) + 1;
  int rstSize = (Nbit + 31) / 32;
  rst.resize(rstSize, 0);
  vec<uint32_t> buffer(rstSize * 3, 0);
  uint32_t *power10 = &buffer[0], *power10new = power10 + rstSize;
  uint32_t *adder = power10new + rstSize;
  power10[0] = 1;
  uint32_t ten = 10;
  for(int i = ssize - 1; ; --i)
  {
    uint32_t a = s[i] - 48;
    mulHint32(adder, rstSize, &a, 1, power10, rstSize);
    addHint32(&rst[0], &rst[0], adder, rstSize);
    if(i <= 0) break;
    mulHint32(power10new, rstSize, power10, rstSize, &ten, 1);
    std::swap(power10new, power10);
  }
}


// [[Rcpp::export]]
IntegerVector testIntString2hugeInt32(std::string &ss)
{
  vec<uint32_t> rst;
  char *s = &ss[0];
  intString2hugeInt32(rst, s, ss.size());


  if(true)
  {
    int i = rst.size() - 1;
    for(; i >= 0 and rst[i] == 0; --i);
    rst.resize(i + 1);
  }
  IntegerVector result(rst.size() * sizeof(uint32_t) * 8);
  int Nbits = sizeof(uint32_t) * 8;
  uint32_t masks[Nbits];
  uint32_t one = 1;
  for(int i = 0; i < Nbits; ++i) masks[i] = one << i;


  for(int i = 0, iend = rst.size(), k = result.size() - 1; i < iend; ++i)
  {
    for(int u = 0, uend = Nbits; u < uend; ++u)
    {
      result[k] = (rst[i] & masks[u]) != 0;
      --k;
    }
  }


  return result;
}


// [[Rcpp::export]]
IntegerVector testMulHintV32(std::string &xv, std::string &yv)
{
  vec<uint32_t> x, y;
  intString2hugeInt32(x, &xv[0], xv.size());
  intString2hugeInt32(y, &yv[0], yv.size());
  if(true)
  {
    int i = x.size() - 1;
    for(; i >= 0 and x[i] == 0; --i);
    x.resize(i + 1);
    i = y.size() - 1;
    for(; i >= 0 and y[i] == 0; --i);
    y.resize(i + 1);
  }


  vec<uint32_t> rst(x.size() + y.size());
  mulHint32(&rst[0], rst.size(), &x[0], x.size(), &y[0], y.size());


  IntegerVector result(rst.size() * sizeof(uint32_t) * 8);
  int Nbits = sizeof(uint32_t) * 8;
  uint32_t masks[Nbits];
  uint32_t one = 1;
  for(int i = 0; i < Nbits; ++i) masks[i] = one << i;


  for(int i = 0, iend = rst.size(), k = result.size() - 1; i < iend; ++i)
  {
    for(int u = 0, uend = Nbits; u < uend; ++u)
    {
      result[k] = (rst[i] & masks[u]) != 0;
      --k;
    }
  }


  return result;
}


/*
void addHint(uint64_t *z, uint64_t *x, uint64_t *y, int size)
{
  z[0] = x[0] + y[0];
  uint64_t carry = z[0] < x[0] or z[0] < y[0];
  for(int i = 1; i < size; ++i)
  {
    z[i] = x[i] + y[i];
    uint64_t carryNew = z[i] < x[i] or z[i] < y[i];
    z[i] += carry;
    carry = carryNew;
  }
}
*/


// [[Rcpp::export]]
List testAddHintVbinInput(IntegerVector xv, IntegerVector yv)
{
  vec<int> x(xv.begin(), xv.end()), y(yv.begin(), yv.end());
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());
  int asize = (x.size() + 63) / 64;
  int bsize = (y.size() + 63) / 64;
  vec<uint64_t> a(asize, 0), b(bsize, 0);
  for(int i = 0, k = 0; i < asize; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      a[i] += ((uint64_t)x[k]) << j;
  }
  for(int i = 0, k = 0; i < bsize; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      b[i] += ((uint64_t)y[k]) << j;
  }


  if(false)
  {
    std::string astr = "";
    for(int i = a.size() - 1; i >= 0; --i)
      astr += std::bitset<64>(a[i]).to_string();
    std::string bstr = "";
    for(int i = b.size() - 1; i >= 0; --i)
      bstr += std::bitset<64>(b[i]).to_string();
    std::cout << "astr = " << astr << std::endl;
    std::cout << "bstr = " << bstr << std::endl;
  }


  vec<uint64_t> ab(a.size(), 0);
  // return List::create();
  addHint(&ab[0], &a[0], &b[0], b.size());
  std::string abstr = "";
  for(int i = ab.size() - 1; i >= 0; --i)
    abstr += std::bitset<64>(ab[i]).to_string();


  vec<int> sum(x.size(), 0);
  for(int i = 0, iend = x.size(); i < iend; ++i)
  {
    sum[i] = x[i] + y[i];
  }
  for(int i = 0, iend = sum.size() - 1; i < iend; ++i)
  {
    sum[i + 1] += sum[i] / 2;
    sum[i] %= 2;
  }


  std::reverse(sum.begin(), sum.end());


  std::string correctRst(sum.size(), '0');
  for(int i = 0, iend = sum.size(); i < iend; ++i)
    correctRst[i] = sum[i] == 1 ? '1' : '0';


  return List::create(Named("rst") = abstr, Named("correctRst") = correctRst);
}


// [[Rcpp::export]]
IntegerVector stringToBinIntVec(std::string &ss)
{
  vec<uint32_t> rst;
  char *s = &ss[0];
  intString2hugeInt32(rst, s, ss.size());


  if(true)
  {
    int i = rst.size() - 1;
    for(; i >= 0 and rst[i] == 0; --i);
    rst.resize(i + 1);
  }


  vec<int> result(rst.size() * sizeof(uint32_t) * 8);
  int Nbits = sizeof(uint32_t) * 8;
  uint32_t masks[Nbits];
  uint32_t one = 1;
  for(int i = 0; i < Nbits; ++i) masks[i] = one << i;


  for(int i = 0, iend = rst.size(), k = 0; i < iend; ++i)
  {
    for(int u = 0, uend = Nbits; u < uend; ++u)
    {
      result[k] = (rst[i] & masks[u]) != 0;
      ++k;
    }
  }


  if(true)
  {
    int i = result.size() - 1;
    for(; i >= 0 and result[i] == 0; --i);
    result.resize(i + 1);
  }


  return IntegerVector(result.begin(), result.end());
}


// [[Rcpp::export]]
int leftTrailingZeros(IntegerVector x)
{
  int i = 0;
  for(int iend = x.size(); i < iend and x[i] != 0; ++i);
  return i;
}


// [[Rcpp::export]]
IntegerVector binize(IntegerVector x)
{
  for(int i = 0, iend = x.size() - 1; i < iend; ++i)
  {
    x[i + 1] += x[i] / 2;
    x[i] %= 2;
  }
  return x;
}


// [[Rcpp::export]]
bool binvecGE(IntegerVector x, IntegerVector y)
{
  int i = x.size() - 1;
  for(; i >= 0 and x[i] == y[i]; --i);
  return i < 0 or x[i] > y[i];
}
// [[Rcpp::export]]
bool binvecLE(IntegerVector x, IntegerVector y)
{
  int i = x.size() - 1;
  for(; i >= 0 and x[i] == y[i]; --i);
  return i < 0 or x[i] < y[i];
}


template<typename T>
void removeHighEndZeros(vec<T> &x)
{
  int i = x.size() - 1;
  for(; i >= 1 and x[i] == 0; --i);
  x.resize(i + 1); // x should be of size 1 at least, if x is nonempty.
}


// x's size can be adjusted.
template<typename T>
void normalizeDecimalIntVec(vec<T> &v)
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
bool nonnegVecLessThan(T *x, int xsize, T *y, int ysize)
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
    for(int i = 0, iend = x.size() - 1; i < iend; ++i)
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
    for(int i = 0, iend = x.size(); i < iend; ++i)
    {
      char *u = x[i];
      int usize = xsizes[i];
      sign[i] = (u[0] != '-') * 2 - 1;
      int k = 0;
      for(; k < usize and u[k] != '.'; ++k);
      Ndeci[i] = std::max(usize - 1 - k, 0);
      maxDecimal = std::max<int> (Ndeci[i], maxDecimal);
    }
    // Rcout << "maxDecimal = " << maxDecimal << "\n";


    vec<int> &zerosToAdd = Ndeci;
    for(int i = 0, iend = zerosToAdd.size(); i < iend; ++i)
    {
      zerosToAdd[i] = maxDecimal - zerosToAdd[i];
      // Rcout << "zerosToAdd[i] = " << zerosToAdd[i] << ", ";
    }
    // Rcout << "\n";
    // return;


    value.resize(zerosToAdd.size());
    for(int i = 0, iend = zerosToAdd.size(); i < iend; ++i)
    {
      char *u = x[i];
      int usize = xsizes[i];
      value[i].reserve(usize + zerosToAdd[i]);
      value[i].assign(zerosToAdd[i], 0);
      for(int k = usize - 1; k >= 0; --k)
      {
        if(u[k] != '-' and u[k] != '.') value[i].push_back(u[k] - 48);
      }
    }


    // for(int i = 0, iend = sign.size(); i < iend; ++i)
    // {
    //   for(int k = 0, kend = value[i].size(); k < kend; ++k)
    //     Rcout << int(value[i][k]);
    //   Rcout << "\n";
    // }
    // Rcout << "\n";
    //
    //
    // for(int i = 0, iend = sign.size(); i < iend; ++i)
    //   Rcout << int(sign[i]) << ", ";
    // Rcout << "\n";
    // Rcout << "\n\n";


    // Find minimum.
    int mini = 0;
    for(int i = 1, iend = sign.size() - 1; i < iend; ++i) // Do not include target.
    {
      vec<int> &minx = value[mini];
      vec<int> &current = value[i];
      bool currentIsLess = valuevecLessThan(
        &current[0], current.size(), &minx[0], minx.size(), sign[i], sign[mini]);
      if(currentIsLess) mini = i;
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
        // Rcout << "\n===================\n";
        for(int k = 0, kend = v.size(); k < kend; ++k)
        {
          // Rcout << v[k] << ", ";
          v[k] = -v[k];
        }
        // Rcout << "\n===================\n";
      }
      // Rcout << "\n++++++++++++++++\n";
      for(int k = 0, kend = std::min(v.size(), vmin.size()); k < kend; ++k)
      {
        v[k] -= vmin[k] * multiple;
        // Rcout << v[k] << ", ";
      }
      // Rcout << "\n++++++++++++++++\n";
      normalizeDecimalIntVec<int> (v);
    }
    value[mini].assign(1, 0);


    targetSum.swap(value.back());
    value.resize(value.size() - 1);


    // valid.resize(value.size());
    // for(int i = 0, iend = value.size(); i < iend; ++i)
    // {
    //   valid[i] = !nonnegVecLessThan(
    //     &targetSum[0], targetSum.size(), &value[i][0], value[i].size());
    // }
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
    std::sort(order.begin(), order.end(), ComparePosiVec(&value[0]));
    rank.resize(order.size());
    for(int i = 0, iend = rank.size(); i < iend; ++i)
      rank[order[i]] = i;
  }


  vec<int> largestSubsetSum;
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


  /*
  vec<vec<int> > stringsToBitIntVec_test(
      StringMatrix &Xval, int whichCol, StringVector &tsum, int len)
  {
    stringsToBitIntVec(Xval, whichCol, tsum, len);
    for(int i = 0, iend = value.size(); i < iend; ++i)
    {
      vec<int> &v = value[i];
      std::reverse(v.begin(), v.end());
    }
    return value;
  }
  */
};


/*
// [[Rcpp::export]]
List AdjustStringVecTest(
    StringMatrix &Xval, StringVector &tsum, int len)
{
  int Ncol = tsum.size();
  vec<vec<vec<int> > > rst(Ncol);
  AdjustStringVec asv;
  vec<vec<int> > targetSums(Ncol);
  for(int i = 0; i < Ncol; ++i)
  {
    rst[i] = asv.stringsToBitIntVec_test(Xval, i, tsum, len);
    targetSums[i] = asv.targetSum;
  }

  return List::create(Named("v") = rst, Named("sum") = targetSums);
}
*/


// Insert at most 64 bits in y into x.
// Return true if y has been exhausted.
void insertBits(uint64_t *x, int &XbeginBit,
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
void insertAllBits(uint64_t *x, vec<vec<uint64_t> > &y,
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


// [[Rcpp::export]]
IntegerVector testInsertBits(List X, IntegerVector NeffectiveBits)
{
  vec<vec<uint64_t> > y(X.size());
  for(int i = 0, iend = y.size(); i < iend; ++i)
  {
    IntegerVector x = X[i];
    y[i].resize(x.size() / 2);
    std::memcpy(&y[i][0], &x[0], sizeof(uint64_t) * y[i].size());
  }
  int totalBits = std::accumulate(
    NeffectiveBits.begin(), NeffectiveBits.end(), 0);
  int finalNint32bit = (totalBits + (sizeof(uint32_t) * 8 - 1)) /
    (sizeof(uint32_t) * 8);
  IntegerVector rst(finalNint32bit, 0);
  insertAllBits((uint64_t*)(&rst[0]), y, &NeffectiveBits[0]);
  return rst;
}




// Return true if the matrix is not empty.
// order is the row order.
bool stringMatTo64bitIntMat(StringMatrix &Xval, StringVector &tsum,
                            int len, int maxCore,
                            vec<uint64_t> &X, vec<uint64_t> &targetSS,
                            vec<uint64_t> &largestSS, vec<int> &order,
                            vec<int> &colOrder, int &totalEffectiveBits)
{
  auto empf = [](std::size_t t){};


  vec<AdjustStringVec> adjv(tsum.size());
  for(int i = 0, iend = adjv.size(); i < iend; ++i)
    adjv[i].stringsToBitIntVecInitialize(Xval, i, tsum, len);


  ParaFor(0, adjv.size(), [&](std::size_t i, std::size_t t)
  {
    adjv[i].stringsToBitIntVec(len);
  }, maxCore, 1, empf, empf);


  // If targetSum of a certain dimension is negative, solution does not exist.
  for(int i = 0, iend = adjv.size(); i < iend; ++i)
  {
    if(adjv[i].targetSum.back() < 0) return false;
  }


  // Compute rank orders and largest subset sums.
  ParaFor(0, adjv.size(), [&](std::size_t i, std::size_t t)
  {
    adjv[i].createRankOrder();
    adjv[i].getLargestSubsetSum(len);
  }, maxCore, 1, empf, empf);


  // Compute rank correlations with others (just inner product)
  vec<std::size_t> innerProd(adjv.size(), 0);
  ParaFor(0, adjv.size(), [&](std::size_t i, std::size_t t)
  {
    auto &x = adjv[i];
    for(int k = 0, kend = adjv.size(); k < kend; ++k)
      innerProd[i] += std::inner_product(
        x.rank.begin(),  x.rank.end(), adjv[k].rank.begin(), 0);
  }, maxCore, 1, empf, empf);


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
  ParaFor(0, adjv.size(), [&](std::size_t i, std::size_t t)
  {
    Nbits[i] = adjv[i].generate64bitUint();
  }, maxCore, 1, empf, empf);


  // # of 64-bit integers in every row.
  int totalBits = std::accumulate(Nbits.begin(), Nbits.end(), 0);
  totalEffectiveBits = totalBits;
  uint64_t dim = (totalBits + sizeof(uint64_t) * 8 - 1) / (sizeof(uint64_t) * 8);
  X.resize(adjv[0].value64bit.size() * dim);


  vec<int> NbitsColOrdered(Nbits.size());
  for(int i = 0, iend = Nbits.size(); i < iend; ++i)
    NbitsColOrdered[i] = Nbits[colOrder[i]];
  Nbits.swap(NbitsColOrdered);


  // Put all the bits in X.
  ParaFor(0, adjv[0].value.size(), [&](std::size_t i, std::size_t t)
  {
    vec<vec<uint64_t> > row(adjv.size());
    for(int j = 0, jend = row.size(); j < jend; ++j)
      row[j].swap(adjv[ colOrder[j] ].value64bit[i]);
    insertAllBits(&X[0] + dim * i, row, &Nbits[0]);
  }, maxCore, 1, empf, empf);


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


  return true;
}




// [[Rcpp::export]]
List stringMatTo64bitIntMatTest(StringMatrix &Xval, StringVector &tsum,
                                int len, int maxCore)
{
  vec<uint64_t> Xbit, targetSS, largestSS;
  vec<int> order, colOrder;
  int totalBits = 0;
  bool success = stringMatTo64bitIntMat(
    Xval, tsum, len, maxCore, Xbit, targetSS,
    largestSS, order, colOrder, totalBits);
  if(!success) return List::create();
  IntegerMatrix XbitR(largestSS.size() * 2, Xbit.size() / largestSS.size());
  IntegerVector targetSSR(targetSS.size() * 2);
  IntegerVector largestSSR(largestSS.size() * 2);
  std::memcpy(&XbitR[0], &Xbit[0], Xbit.size() * sizeof(uint64_t));
  std::memcpy(&targetSSR[0], &targetSS[0], targetSS.size() * sizeof(uint64_t));
  std::memcpy(&largestSSR[0], &largestSS[0], largestSS.size() * sizeof(uint64_t));
  return List::create(
    Named("NeffectiveBitInRow") = totalBits,
    Named("order") = IntegerVector(order.begin(), order.end()) + 1,
    Named("Xbit") = transpose(XbitR), Named("targetBit") = targetSSR,
    Named("largestSSbit") = largestSSR,
    Named("colOrder") = IntegerVector(colOrder.begin(), colOrder.end()) + 1);
}




































/*
// Every 64bit integer has a carry bit.
void Hint32to64(vec<uint64_t> &rst, uint32_t *x, int size)
{
  int Nint64 = (size * 32 + 62) / 63;
  rst.resize(Nint64, 0);
  vec<uint64_t> y(Nint64, 0);
  std::memcpy(&y[0], x, 32 * size);
  for(int i = 0; i < Nint64; ++i)
  {
    int s = i * 63;
    int whichInt = s / 64, whichBit = s % 64;
    uint64_t lower = y[whichInt] >> whichBit;
    whichInt = (s + 62) / 64, whichBit = (s + 62) % 64; // last bit position
    uint64_t upper = y[whichInt] << (63 - whichBit);
    rst[i] = lower + upper;
  }
}


constexpr uint64_t ENDBIT = ((uint64_t)1) << (sizeof(uint64_t) * 8 - 1);
constexpr uint64_t _ENDBIT = ~ENDBIT;
constexpr uint64_t Nfull = sizeof(uint64_t) * 8;
constexpr uint64_t Nfull_1 = Nfull - 1;


// z = x + y, x or/and y can be z.
void addHint(uint64_t *z, uint64_t *x, uint64_t *y, int size)
{
  z[0] = x[0] + y[0];
  for(int i = 1; i < size; ++i)
  {
    uint64_t carry = z[i - 1] >> Nfull_1;
    z[i - 1] &= _ENDBIT;
    z[i] = x[i] + y[i] + carry;
  }
}


// [[Rcpp::export]]
List testAddHintV(IntegerVector xv, IntegerVector yv)
{
  vec<int> x(xv.begin(), xv.end()), y(yv.begin(), yv.end());
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());
  int size = x.size() / 64;
  vec<uint64_t> a(size, 0), b(size, 0);
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      a[i] += ((uint64_t)x[k]) << j;
  }
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      b[i] += ((uint64_t)y[k]) << j;
  }


  if(false)
  {
    std::string astr = "";
    for(int i = a.size() - 1; i >= 0; --i)
      astr += std::bitset<64>(a[i]).to_string();
    std::string bstr = "";
    for(int i = b.size() - 1; i >= 0; --i)
      bstr += std::bitset<64>(b[i]).to_string();
    std::cout << "astr = " << astr << std::endl;
    std::cout << "bstr = " << bstr << std::endl;
  }


  vec<uint64_t> ab(a.size(), 0);
  addHint(&ab[0], &a[0], &b[0], a.size());
  // mulHint(&ab[0], ab.size(), &a[0], a.size(), &b[0], b.size());
  std::string abstr = "";
  for(int i = ab.size() - 1; i >= 0; --i)
    abstr += std::bitset<64>(ab[i]).to_string();


  vec<int> xnew(x.size(), 0), ynew(y.size(), 0);
  // Skip every 64th bit.
  for(int i = 0, iend = x.size() / 64, *xnewk = &xnew[0]; i < iend; ++i)
  {
    int *lb = &x[0] + i * 64;
    int *ub = lb + 63;
    std::copy(lb, ub, xnewk);
    xnewk += 63;
  }
  for(int i = 0, iend = y.size() / 64, *ynewk = &ynew[0]; i < iend; ++i)
  {
    int *lb = &y[0] + i * 64;
    int *ub = lb + 63;
    std::copy(lb, ub, ynewk);
    ynewk += 63;
  }
  x.swap(xnew); y.swap(ynew);


  vec<int> added(x.size(), 0);
  for(int i = 0, iend = x.size(); i < iend; ++i) added[i] = x[i] + y[i];


  for(int i = 0, iend = added.size() - 1; i < iend; ++i)
  {
    added[i + 1] += added[i] / 2;
    added[i] %= 2;
  }
  vec<int> addednew(added.size() * 2, 0);
  for(int i = 0, iend = (added.size() + 62) / 63; i < iend; ++i)
  {
    int *lb = &added[0] + i * 63;
    int *ub = std::min(lb + 63, &*added.end());
    std::copy(lb, ub, addednew.begin() + i * 64);
  }
  addednew.resize(added.size());
  added.swap(addednew);
  added.back() = added[added.size()];


  std::reverse(added.begin(), added.end());


  std::string correctRst(added.size(), '0');
  for(int i = 0, iend = added.size(); i < iend; ++i)
    correctRst[i] = added[i] == 1 ? '1' : '0';


  return List::create(Named("rst") = abstr, Named("correctRst") = correctRst);
}


// z = x - y, x or/and y can be z.
void subHint(uint64_t *z, uint64_t *x, uint64_t *y, int size)
{
  z[0] = x[0] - y[0];
  for(int i = 1; i < size; ++i)
  {
    uint64_t carry = z[i - 1] >> Nfull_1;
    z[i - 1] &= _ENDBIT;
    z[i] = x[i] - carry - y[i];
  }
}


// [[Rcpp::export]]
List testSubHintV(IntegerVector xv, IntegerVector yv)
{
  vec<int> x(xv.begin(), xv.end()), y(yv.begin(), yv.end());
  std::reverse(x.begin(), x.end());
  std::reverse(y.begin(), y.end());
  int size = x.size() / 64;
  vec<uint64_t> a(size, 0), b(size, 0);
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      a[i] += ((uint64_t)x[k]) << j;
  }
  for(int i = 0, k = 0; i < size; ++i)
  {
    for(int j = 0; j < 64; ++k, ++j)
      b[i] += ((uint64_t)y[k]) << j;
  }


  if(false)
  {
    std::string astr = "";
    for(int i = a.size() - 1; i >= 0; --i)
      astr += std::bitset<64>(a[i]).to_string();
    std::string bstr = "";
    for(int i = b.size() - 1; i >= 0; --i)
      bstr += std::bitset<64>(b[i]).to_string();
    std::cout << "astr = " << astr << std::endl;
    std::cout << "bstr = " << bstr << std::endl;
  }


  vec<uint64_t> ab(a.size(), 0);
  subHint(&ab[0], &a[0], &b[0], a.size());
  // mulHint(&ab[0], ab.size(), &a[0], a.size(), &b[0], b.size());
  std::string abstr = "";
  for(int i = ab.size() - 1; i >= 0; --i)
    abstr += std::bitset<64>(ab[i]).to_string();


  vec<int> xnew(x.size(), 0), ynew(y.size(), 0);
  // Skip every 64th bit.
  for(int i = 0, iend = x.size() / 64, *xnewk = &xnew[0]; i < iend; ++i)
  {
    int *lb = &x[0] + i * 64;
    int *ub = lb + 63;
    std::copy(lb, ub, xnewk);
    xnewk += 63;
  }
  for(int i = 0, iend = y.size() / 64, *ynewk = &ynew[0]; i < iend; ++i)
  {
    int *lb = &y[0] + i * 64;
    int *ub = lb + 63;
    std::copy(lb, ub, ynewk);
    ynewk += 63;
  }
  x.swap(xnew); y.swap(ynew);


  vec<int> subed(x.size(), 0);
  for(int i = 0, iend = x.size(); i < iend; ++i) subed[i] = x[i] - y[i];


  for(int i = 0, iend = subed.size() - 1; i < iend; ++i)
  {
    // subed[i + 1] += subed[i] / 2;
    subed[i + 1] += (subed[i] - 1) / 2;
    subed[i] = std::abs(subed[i] % 2);
  }
  vec<int> subednew(subed.size() * 2, 0);
  for(int i = 0, iend = (subed.size() + 62) / 63; i < iend; ++i)
  {
    int *lb = &subed[0] + i * 63;
    int *ub = std::min(lb + 63, &*subed.end());
    std::copy(lb, ub, subednew.begin() + i * 64);
  }
  subednew.resize(subed.size());
  subed.swap(subednew);
  subed.back() = subed[subed.size()];


  std::reverse(subed.begin(), subed.end());


  std::string correctRst(subed.size(), '0');
  for(int i = 0, iend = subed.size(); i < iend; ++i)
    correctRst[i] = subed[i] == 1 ? '1' : '0';


  return List::create(Named("rst") = abstr, Named("correctRst") = correctRst);
}
*/







































