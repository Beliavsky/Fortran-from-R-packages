// -*- C -*-
#include "randompack_config.h"

#include <stdlib.h>
#include <stdbool.h>

// let:
//   A = [a0..a8]
//   B = [b0..b8]
//   R = [r0..r17]
// where ai, bi, ri are 64bit integers. Also let
// hi:lo denote a 128 bit integer
// c be 64 bits
// cy be a carry bit
// cout = adc(x, y, cin) be an add-with-carry function that adds y to x

// And the algorithm:
// R = 0
// for i = 0..8:
//   c = 0
//   cy = 0
//   for j = 0..8:
//     hi:lo = b[i]*a[j]
//     hi:lo += c        // 128-bit add
//     c = hi
//     cy = adc(r[i+j], lo, cy)
//   cy = adc(r[i+9], c, cy)

#if defined(_MSC_VER)
ALWAYS_INLINE uint64_t addc64(uint64_t a, uint64_t b, uint64_t cin, uint64_t *cout) {
  uint64_t result;
  *cout = _addcarry_u64((unsigned char)cin, a, b, &result);
  return result;
}
ALWAYS_INLINE uint64_t subb64(uint64_t a, uint64_t b, uint64_t bin, uint64_t *bout) {
  uint64_t result;
  *bout = _subborrow_u64((unsigned char)bin, a, b, &result);
  return result;
}
#else
ALWAYS_INLINE uint64_t addc64(uint64_t a, uint64_t b, uint64_t cin, uint64_t *cout) {
  uint64_t s;
  bool c0 = __builtin_add_overflow(a, b, &s);
  bool c1 = __builtin_add_overflow(s, cin, &s);
  *cout = (uint64_t)(c0 || c1);
  return s;
}
ALWAYS_INLINE uint64_t subb64(uint64_t a, uint64_t b, uint64_t bin, uint64_t *bout) {
  uint64_t d;
  bool b0 = __builtin_sub_overflow(a, b, &d);
  bool b1 = __builtin_sub_overflow(d, bin, &d);
  *bout = (uint64_t)(b0 || b1);
  return d;
}
#endif

static const uint64_t ranlux_a[9] = {
  0xed7faa90747aaad9ULL,
  0x4cec2c78af55c101ULL,
  0xe64dcb31c48228ecULL,
  0x6d8a15a13bee7cb0ULL,
  0x20b2ca60cb78c509ULL,
  0x256c3d3c662ea36cULL,
  0xff74e54107684ed2ULL,
  0x492edfcc0cc8e753ULL,
  0xb48c187cf5b22097ULL,
};

ALWAYS_INLINE void mul9x9(uint64_t *z, const uint64_t *restrict x,
  const uint64_t *restrict y) {
  uint64_t lo, hi, c, cy;
  for (int i=0; i<18; i++) z[i] = 0;
#if defined(__clang__) || (defined(__GNUC__) && !defined(__NVCOMPILER))
#pragma GCC unroll 9
#endif
  for (int i=0; i<9; i++) {
    c = 0;
    cy = 0;
    for (int j=0; j<9; j++) {
      MUL64_WIDE(x[i], y[j], hi, lo);
      add128(lo, hi, c, 0, &lo, &hi);
      c = hi;
      z[i+j] = addc64(z[i+j], lo, cy, &cy);
    }
    z[i+9] = c + cy;
  }
}

// Remainder algorithm:
// [let m = 2^576 - 2^240 + 1]
// t0 = z[0..8]
// t1 = z[9..17]
// t2 = upper 240 bits of t1 (3.75 words)
// t3 = lower 336 bits of t1 (5.25 words)
// r = t0 - t1 - t2 + (t3 + t2)*2^240
// c = floor(r/2^576)
// x = r - c*m
// [it is easy to see that 0 <= x < m]

ALWAYS_INLINE int add1_cy(uint64_t *x, uint64_t c, int n) {
  uint64_t carry = 0;
  x[0] = addc64(x[0], c, 0, &carry);
  for (int i = 1; i < n; i++) {
    x[i] = addc64(x[i], 0, carry, &carry);
  }
  return (int)carry;
}

ALWAYS_INLINE int sub1_bw(uint64_t *x, uint64_t c, int n) {
  uint64_t borrow = 0;
  x[0] = subb64(x[0], c, 0, &borrow);
  for (int i = 1; i < n; i++) x[i] = subb64(x[i], 0, borrow, &borrow);
  return (int)borrow;
}

ALWAYS_INLINE int add_cy9(uint64_t *r, const uint64_t *restrict y) {
  uint64_t carry = 0;
  for (int i = 0; i < 9; i++) r[i] = addc64(r[i], y[i], carry, &carry);
  return (int)carry;
}

ALWAYS_INLINE int add_cy4(uint64_t *r, const uint64_t y[9]) {
  uint64_t carry = 0;
  for (int i = 0; i < 5; i++) r[i] = addc64(r[i], y[i], carry, &carry);
  if (!carry) return 0;
  for (int i = 5; i < 9; i++) r[i] = addc64(r[i], 0, carry, &carry);
  return (int)carry;
}

ALWAYS_INLINE int sub_bw9(uint64_t *r, const uint64_t *restrict y) {
  uint64_t borrow = 0;
  for (int i = 0; i < 9; i++) r[i] = subb64(r[i], y[i], borrow, &borrow);
  return (int)borrow;
}

ALWAYS_INLINE int sub_bw4(uint64_t *r, const uint64_t y[9]) {
  uint64_t borrow = 0;
  for (int i = 0; i < 5; i++) r[i] = subb64(r[i], y[i], borrow, &borrow);
  if (!borrow) return 0;
  for (int i = 5; i < 9; i++) r[i] = subb64(r[i], 0, borrow, &borrow);
  return (int)borrow;
}

ALWAYS_INLINE void t2get(uint64_t t2[4], uint64_t t1[9]) {
  t2[0] = (t1[5] >> 16) | (t1[6] << 48);
  t2[1] = (t1[6] >> 16) | (t1[7] << 48);
  t2[2] = (t1[7] >> 16) | (t1[8] << 48);
  t2[3] =  t1[8] >> 16;
}

ALWAYS_INLINE int shly240(uint64_t y[9]) {
  int co = (int)((y[5] >> 16) & 1);
  y[8] = (y[5] << 48) | (y[4] >> 16);
  y[7] = (y[4] << 48) | (y[3] >> 16);
  y[6] = (y[3] << 48) | (y[2] >> 16);
  y[5] = (y[2] << 48) | (y[1] >> 16);
  y[4] = (y[1] << 48) | (y[0] >> 16);
  y[3] =  y[0] << 48;
  y[2] = 0;
  y[1] = 0;
  y[0] = 0;
  return co;
}

ALWAYS_INLINE void corr_r_minus_cm(uint64_t r[9], int c) {
  uint64_t
    c0 = (uint64_t)abs(c),
    c3 = ((uint64_t)abs(c)) << 48;
  if (c > 0) {
    sub1_bw(r, c0, 9);
    add1_cy(r + 3, c3, 6);
  }
  else if (c < 0) {
    add1_cy(r, c0, 9);
    sub1_bw(r + 3, c3, 6);
  }
}

ALWAYS_INLINE void mod9x9(uint64_t x[9], uint64_t z[18]) {
  uint64_t t1[9], t2[9] = {0}, t3[9];
  int c;
  copy64(t1, z + 9, 9);
  t2get(t2, t1);
  copy64(t3, z + 9, 6); t3[5] &= 0xFFFF;
  copy64(x, z, 9);
  c = 0;
  c -= sub_bw9(x, t1);
  c -= sub_bw4(x, t2);
  add_cy4(t3, t2);
  shly240(t3);
  c += add_cy9(x, t3);
  corr_r_minus_cm(x, c);
}

static inline ALIGN64 void fill_ranluxpp(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  uint64_t s[9];
  uint64_t z[18];
  copy64(s, state->u64, 9);
  size_t i = 0;
  while (i < len) {
    mul9x9(z, s, ranlux_a);
    mod9x9(s, z);
    for (int j = 0; j < 9 && i < len; j++, i++)
      out[i] = s[j];
  }
  copy64(state->u64, s, 9);
}
