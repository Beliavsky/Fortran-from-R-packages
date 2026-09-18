// -*- C -*-
// Code copied from NumPy: github.com/numpy/numpy/tree/f6440be7b8,
//   numpy/random/src/distributions/distributions.c.
// and adapted to randompack. No change to computations.
// Copyright (c) NumPy Developers.
//
// Strips are numbered from top to bottom: 1, 2,...,255, 0
// Bottom right corner of strip i is at (xi, fi).
// Accept if x <= x{i-1}.
// Draw y if x{i-1} < x <= x{i}

#define DBL(x) ((x)*0x1.0p-52)

#include <stdint.h>
#include <math.h>
#include <string.h>
#include "printX.h"
#include "randompack.h"
#include "ziggurat_const_double.h"
#include "zig_D.h"
#include "assert.h"

//=============================== NORMAL DISTRIBUTION ==============================

ALWAYS_INLINE double norm_decode(uint64_t r, uint8_t idx, uint64_t *ki,
  uint64_t *rabs, bool *sign) {
  *ki = ki_double[idx];
  *rabs = (r >> 9) & 0x000fffffffffffff;
  *sign = ((r >> 8) & 1) != 0;
  int64_t s = (int64_t)*rabs;
#if defined(_MSC_VER) && !defined(__clang__)
  int64_t m = -(int64_t)*sign;
  s = (s ^ m) - m;
#else
  if (*sign) s = -s;
#endif
  return (double)s * wi_double[idx];
}

NEVER_INLINE bool below_normpdf(uint8_t idx, double x, double y) {
  double f0 = fi_double[idx-1];
  double f1 = fi_double[idx];
  return f1 + y*(f0 - f1) < openlibm_exp(-x*x/2);
}

NEVER_INLINE double norm_tail_draw(randompack_rng *rng, bool sign, double y, int *w) {
  double xx, yy, x;
  for (;;) {
    /* Switch to 1.0 - U to avoid log(0.0), see GH 13361 */
    x = draw_u01(rng, w);
    xx = -ziggurat_nor_inv_r * openlibm_log1p(-x);
    yy = -openlibm_log1p(-y);
    if (LIKELY(yy + yy > xx * xx))
      return sign ? -(ziggurat_nor_r + xx) : ziggurat_nor_r + xx;
    y = draw_u01(rng, w);
  }
}

NEVER_INLINE double norm_unlikely(randompack_rng *rng, int *w, double x, uint8_t idx,
  uint64_t rabs, uint64_t ki, bool sign) {
  uint64_t ru, tmp, L, R, Y, YL;
  double u = 0;
  const int inflection_idx = 52;
  for (;;) {
    Y = draw_u64(rng, w);
    if (UNLIKELY(idx == 0)) {
      u = U64_TO_DOUBLE(Y);
      x = norm_tail_draw(rng, sign, u, w);
      return x;
    }
    // The following block of 16 lines (until the blk-end marker) has been checked to be
    // correct by comparing draws of ~10^14 values with it commented out and commented in.
    L = (1ULL << 52) - ki;
    R = (1ULL << 52) - rabs;
    MUL64_WIDE(Y, L, YL, tmp);  //
    (void) tmp;
    if (LIKELY(idx > inflection_idx)) {
      if (YL > R) goto next_try;
      if (YL + Dnorm[idx] < R) return x;
    }
    else if (LIKELY(idx < inflection_idx)) {
      if (YL < R) return x;
      if (YL > R + Dnorm[idx]) goto next_try;
    }
    else {
      if (YL > R + Dzig52) goto next_try;
      if (YL + Dnorm[idx] < R) return x;
    } // blk-end
    u = U64_TO_DOUBLE(Y);
    if (below_normpdf(idx, x, u)) return x;
  next_try:
    ru = draw_u64(rng, w);
    idx = ru & 0xff;
    x = norm_decode(ru, idx, &ki, &rabs, &sign);
    if (LIKELY(rabs < ki)) return x;
  }
}

ALWAYS_INLINE double draw_norm(randompack_rng *rng, int *w) {
  uint64_t rabs, ki, r;
  double x;
  uint8_t idx;
  bool sign;
  r = draw_u64(rng, w);
  idx = r & 0xff;
  x = norm_decode(r, idx, &ki, &rabs, &sign);
  if (LIKELY(rabs < ki)) return x;
  x = norm_unlikely(rng, w, x, idx, rabs, ki, sign);
  return x;
}

NEVER_INLINE void rand_norm(double *buf, size_t len, randompack_rng *rng) {
  align64(rng);
  draw_raw_copy((uint8_t *)buf, len*sizeof(uint64_t), rng);
  int w = enter_u64_mode(rng);
  uint8_t idx;
  uint64_t r, rabs, ki;
  double x;
  bool sign;
  for (size_t i = len; i-- > 0;) {
    uint8_t *bufpi = (uint8_t*)(buf + i);
    memcpy(&r, bufpi, sizeof(r));
    idx = r & 0xff;
    x = norm_decode(r, idx, &ki, &rabs, &sign);
    buf[i] = x;
    buf[i] = LIKELY(rabs < ki) ? x : norm_unlikely(rng, &w, x, idx, rabs, ki, sign);
  }
  exit_u64_mode(rng, w);
}

//========================= exponential distribution ===========================

ALWAYS_INLINE bool below_exppdf(uint8_t idx, double x, double y) {
  // Use openlibm_exp to ensure cross-platform compatibility (with exp we
  // the result of the test might depend on the platform)
  double f0, f1;
  f0 = fe_double[idx-1];
  f1 = fe_double[idx];
  return (f1 + y*(f0 - f1) < openlibm_exp(-x)); // below curve?
}

ALWAYS_INLINE double exp_taildraw(randompack_rng *rng, int *w) {
  double x = ziggurat_exp_r - openlibm_log1p(-draw_u01(rng, w));
  return x;
}

NEVER_INLINE double exp_unlikely(randompack_rng *rng, int *w, double x, uint8_t idx,
  uint64_t r) {
  double y;
  uint64_t Y, YL, R, L, tmp, ke = ke_double[idx];
  while (true) {
    if (UNLIKELY(idx == 0)) return exp_taildraw(rng, w);
    Y = draw_u64(rng, w);
    L = (1ULL << 53) - ke;
    R = (1ULL << 53) - r;
    MUL64_WIDE(Y, L, YL, tmp);
    (void)tmp;
    if (LIKELY(YL > R)) goto next_try;         // above secant
    if (YL + Dexp[idx] < R) return x;  // below pdf by the gap bound
    y = U64_TO_DOUBLE(Y);
    if (below_exppdf(idx, x, y)) return x; // below curve
  next_try:
    r = draw_u64(rng, w);
    idx = r & 0xff;
    r = (r >> 8) & 0x1FFFFFFFFFFFFF;
    ke = ke_double[idx];
    x = r*we_double[idx];
    if (LIKELY(r < ke)) return x;
  }
}

NEVER_INLINE void rand_exp(double *buf, size_t len, randompack_rng *rng) {
  align64(rng);
  draw_raw_copy((uint8_t *)buf, len*sizeof(uint64_t), rng);
  int w = enter_u64_mode(rng);
  uint8_t idx;
  uint64_t rraw, r;
  double x;
  for (double *p = buf + len; p-- > buf; ) {
    memcpy(&rraw, p, sizeof(rraw));
    idx = rraw & 0xff;
    r = (rraw >> 8) & 0x1FFFFFFFFFFFFF;
    x = r*we_double[idx];
    *p = LIKELY(r < ke_double[idx]) ? x : exp_unlikely(rng, &w, x, idx, r);
  }
  exit_u64_mode(rng, w);
}
