#pragma once
#include <cstdint>


#if __has_include("gmp.h")
#define hasgmp
#include <gmp.h>
#endif


#ifndef vec
#define vec std::vector
#endif


// z can equal x or y.
inline void addHint32(uint32_t *z, uint32_t *x, uint32_t *y, int size)
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
inline void mulHint32(uint32_t *z, int zsize,
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
      // s is the carry;
      for(++k; s != 0; ++k, s >>= 32)
      {
        s += z[k];
        z[k] = uint32_t(s);
      }
    }
  }
}


// The most immediate 8-byte aligned address.
inline void *N8BAA(void *x)
{
  constexpr const uint64_t a = ((uint64_t(0) - 1) >> 3) << 3;
  return reinterpret_cast<void*>((reinterpret_cast<std::uintptr_t>(x) + 7) & a);
}


template<typename T>
inline void subset(vec<T> &X, bool *B, int worthItems)
{
  vec<T> Y(worthItems);
  int j = 0;
  for (int i = 0, iend = X.size(); i < iend; ++i)
  {
    if(B[i])
    {
      std::swap(X[i], Y[j]);
      ++j;
    }
  }
  std::swap(X, Y);
}


// DO NOT ATTEMPT TO INSERT CARRYING BITS ANY MORE!
// THOUGHT AND IMPLEMENTED IT AGAIN, BUT IT TURNS OUT SUBTRACTION'S ABSOLUTE
// VALUE DEMANDS MORE BITWISE OPERATIONS!
// Add huge integers: z[] = x[] + y[]
template<typename ing>
inline void myAddAlgo(uint64_t *z, uint64_t *x, uint64_t *y, ing d)
{
  uint64_t carry = 0;
  for (ing i = 0; i < d; ++i)
  {
    uint64_t zi = x[i] + y[i];
    uint64_t carryNew = zi < x[i];
    z[i] = zi + carry;
    carry = carryNew or z[i] < zi;
  }
}


template<typename ing>
inline void addHint(uint64_t *z, uint64_t *x, uint64_t *y, ing d)
{
#ifdef hasgmp
  // Rcpp::Rcout << "hey it has gmp\n";
  if (sizeof(uintptr_t) == 8) mpn_add_n ( // If the platform is 64-bit:
    (mp_limb_t*)z, (mp_limb_t*)x, (mp_limb_t*)y, d);
  else myAddAlgo<ing> (z, x, y, d);
#else
  myAddAlgo<ing> (z, x, y, d);
#endif
}


template<typename ing>
inline void mySubAlgo(uint64_t *x, uint64_t *z, uint64_t *y, ing d)
{
  uint64_t carry = 0;
  for (ing i = 0; i < d; ++i)
  {
    uint64_t xi = z[i] - y[i];
    uint64_t carryNew = z[i] < xi;
    x[i] = xi - carry;
    carry = carryNew or xi < x[i];
  }
}


// Subtract huge integers: x[] = z[] - y[]
template<typename ing>
inline void subHint(uint64_t *x, uint64_t *z, uint64_t *y, ing d)
{
#ifdef hasgmp
  if (sizeof(uintptr_t) == 8) mpn_sub_n (
    (mp_limb_t*)x, (mp_limb_t*)z, (mp_limb_t*)y, d);
  else mySubAlgo<ing> (x, z, y, d);
#else
  mySubAlgo<ing> (x, z, y, d);
#endif
}


// z = x + y; x, y, z can point to the same address.
template<typename ing>
inline void mvalPlus(uint64_t *z, uint64_t *x, uint64_t *y, ing d)
{
  if (d == 1) { *z = *x + *y; return; }
  addHint<ing> (z, x, y, d);
}


template<typename ing>
inline void iterSum(uint64_t *rst, uint64_t **V, ing *ele, ing eleSize, ing d)
{
  for (ing i = 0; i < eleSize; ++i)
  {
    uint64_t *val = V[ele[i]];
    mvalPlus<ing> (rst, rst, val, d);
  }
}


template<typename ing>
inline void mvalMinus(uint64_t *x, uint64_t *z, uint64_t *y, ing d)
{
  if (d == 1) { *x = *z - *y; return; }
  subHint<ing> (x, z, y, d);
}


// a = x + y - z
template<typename ing>
inline void  mvalPlusMinus(uint64_t *a, uint64_t *x, uint64_t *y,
                           uint64_t *z, ing d)
{
  if (d == 1) { *a = *x + *y - *z; return; }
  addHint<ing> (a, x, y, d);
  subHint<ing> (a, a, z, d);
}


// a = x - y + z
template<typename ing>
inline void mvalMinusPlus(uint64_t *a, uint64_t *x, uint64_t *y,
                          uint64_t *z, ing d)
{
  mvalPlusMinus<ing> (a, x, z, y, d);
}


template<typename ing>
inline bool equal(uint64_t *x, uint64_t *y, ing d)
{
  for(--d; d >= 0 and x[d] == y[d]; --d);
  return d < 0;
}


template<typename ing>
inline bool lessEqual(uint64_t *x, uint64_t *y, ing d)
{
  for (--d; d >= 0 and x[d] == y[d]; --d);
  return d < 0 or x[d] < y[d];
}


template<typename ing>
inline bool greaterEqual(uint64_t *x, uint64_t *y, ing d)
{
  return lessEqual(y, x, d);
}


template<typename ing>
inline uint64_t **mvalLowerBoundLr (
    uint64_t **begin, uint64_t **end, uint64_t *val, ing d)
{
  for (; begin < end; ++begin)
  {
    if ( greaterEqual<ing> (*begin, val, d) ) break;
  }
  return begin;
}


template<typename ing>
inline uint64_t **mvalUpperBoundLr (
    uint64_t **begin, uint64_t **end, uint64_t *val, ing d)
{
  for(--end; end >= begin; --end)
  {
    if( lessEqual<ing> (*end, val, d) ) break;
  }
  return end;
}


template<typename ing>
inline uint64_t isNegative(uint64_t *x, ing d)
{
  return x[d - 1] >> (64 - 1);
}













