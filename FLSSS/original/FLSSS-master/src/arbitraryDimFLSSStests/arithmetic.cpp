#include <Rcpp.h>
using namespace Rcpp;
#include "../arbitraryDimFLSSS/arithmetic.hpp"


#ifndef vec
#define vec std::vector
#endif


template<typename T>
void normalizeBinaryIntVec(vec<T> &v)
{
  for(int k = 0, kend = v.size() - 1; k < kend; ++k)
  {
    if(v[k] >= 0)
    {
      v[k + 1] += v[k] / 2;
      v[k] %= 2;
    }
    else
    {
      T subtractor = (-v[k] + 1) / 2;
      v[k + 1] -= subtractor;
      v[k] += subtractor * 2;
    }
  }
  while(v.back() > 1)
  {
    int vk = v.back();
    v.push_back(vk / 2);
    v[v.size() - 2] = vk % 2;
  }
}


void pushInBit(vec<uint64_t> &x, int *b, int bsize)
{
  int size = (bsize + 63) / 64;
  x.assign(size, 0);
  for(int i = 0, k = 0; i < size; ++i)
  {
    uint64_t sh = 0;
    constexpr uint64_t one = 1;
    for(int kend = std::min<int> ((i + 1) * 64, bsize);
        k < kend; ++k, ++sh)
      x[i] += b[k] * (one << sh);
  }
}


// [[Rcpp::export]]
List addHintTest(IntegerVector xbin, IntegerVector ybin)
{
  vec<int> zbin(xbin.size() + 1);
  zbin.resize(xbin.size());
  for(int i = 0, iend = xbin.size(); i < iend; ++i)
    zbin[i] = xbin[i] + ybin[i];
  normalizeBinaryIntVec(zbin);
  vec<uint64_t> zbin64int;
  pushInBit(zbin64int, &zbin[0], zbin.size());
  IntegerVector correctRst(zbin64int.size() * 2, 0);
  std::memcpy(&correctRst[0], &zbin64int[0], zbin64int.size() * 8);


  vec<uint64_t> x64b, y64b, z64b;
  pushInBit(x64b, &xbin[0], xbin.size());
  pushInBit(y64b, &ybin[0], ybin.size());
  z64b.resize(x64b.size() + 1);
  addHint<int> (&z64b[0], &x64b[0], &y64b[0], x64b.size());
  while(z64b.back() == 0) z64b.pop_back();
  IntegerVector rst(z64b.size() * 2, 0);
  std::memcpy(&rst[0], &z64b[0], z64b.size() * 8);


  return List::create(Named("correctRst") = correctRst, Named("rst") = rst);
}


// [[Rcpp::export]]
List subHintTest(IntegerVector xbin, IntegerVector ybin)
{
  vec<int> zbin(xbin.size() + 1);
  zbin.resize(xbin.size());
  for(int i = 0, iend = xbin.size(); i < iend; ++i)
    zbin[i] = xbin[i] - ybin[i];
  normalizeBinaryIntVec(zbin);
  vec<uint64_t> zbin64int;
  pushInBit(zbin64int, &zbin[0], zbin.size());
  while(zbin64int.back() == 0) zbin64int.pop_back();
  IntegerVector correctRst(zbin64int.size() * 2, 0);
  std::memcpy(&correctRst[0], &zbin64int[0], zbin64int.size() * 8);


  vec<uint64_t> x64b, y64b, z64b;
  pushInBit(x64b, &xbin[0], xbin.size());
  pushInBit(y64b, &ybin[0], ybin.size());
  z64b.resize(x64b.size() + 1);
  subHint<int> (&z64b[0], &x64b[0], &y64b[0], x64b.size());
  while(z64b.back() == 0) z64b.pop_back();
  IntegerVector rst(z64b.size() * 2, 0);
  std::memcpy(&rst[0], &z64b[0], z64b.size() * 8);


  return List::create(Named("correctRst") = correctRst, Named("rst") = rst);
}









































