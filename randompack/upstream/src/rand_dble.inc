// -*- C -*-

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#if HAVE_NEON
#include <arm_neon.h>
#endif

// Convert uint64 RNG output to doubles in [0,1) using the 52-bit mantissa method.
// Each value uses the top 52 bits of the integer, forces the IEEE-754 exponent
// to 1023 (giving a number in [1,2)), then subtracts 1. Results are exactly
// k*2^-52 for k=0..2^52−1 and are bit-identical across all implementations.
//
// When available, this uses SIMD: AArch64 NEON uses two 2-lane chains
// (ILP-friendly) per 4 outputs. Scalar tail handles the last 1-3 values when
// len is not a multiple of 4.
ALWAYS_INLINE void rand_dble_4x52(double x[], size_t len, uint64_t u64[]) {
  size_t i = 0;
#if HAVE_NEON
  uint64x2_t expo = vdupq_n_u64(0x3ff0000000000000ULL);
  float64x2_t one = vdupq_n_f64(1);
  for (; i + 3 < len; i += 4) {
    uint64x2_t r0 = vld1q_u64(u64 + i);
    uint64x2_t r1 = vld1q_u64(u64 + i + 2);
    r0 = vshrq_n_u64(r0, 12);
    r1 = vshrq_n_u64(r1, 12);
    r0 = vorrq_u64(r0, expo);
    r1 = vorrq_u64(r1, expo);
    float64x2_t d0 = vreinterpretq_f64_u64(r0);
    float64x2_t d1 = vreinterpretq_f64_u64(r1);
    d0 = vsubq_f64(d0, one);
    d1 = vsubq_f64(d1, one);
    vst1q_f64(x + i, d0);
    vst1q_f64(x + i + 2, d1);
  }
#else
  for (; i + 3 < len; i += 4) {
    uint64_t b0 = (u64[i] >> 12) | 0x3ff0000000000000ULL;
    uint64_t b1 = (u64[i + 1] >> 12) | 0x3ff0000000000000ULL;
    uint64_t b2 = (u64[i + 2] >> 12) | 0x3ff0000000000000ULL;
    uint64_t b3 = (u64[i + 3] >> 12) | 0x3ff0000000000000ULL;
    double d0, d1, d2, d3;
    memcpy(&d0, &b0, sizeof(d0));
    memcpy(&d1, &b1, sizeof(d1));
    memcpy(&d2, &b2, sizeof(d2));
    memcpy(&d3, &b3, sizeof(d3));
    x[i] = d0 - 1;
    x[i + 1] = d1 - 1;
    x[i + 2] = d2 - 1;
    x[i + 3] = d3 - 1;
  }
#endif
  for (; i < len; i++) {
    uint64_t bits = (u64[i] >> 12) | 0x3ff0000000000000ULL;
    double d;
    memcpy(&d, &bits, sizeof(d));
    x[i] = d - 1;
  }
}

ALWAYS_INLINE void rand_unif_4x52(double x[], size_t len, uint64_t u64[],
  double a, double b) {
  double shift = a - b;
  size_t i = 0;
#if HAVE_NEON
  uint64x2_t expo = vdupq_n_u64(0x3ff0000000000000ULL);
  float64x2_t s = vdupq_n_f64(b);
  float64x2_t c = vdupq_n_f64(shift);
  for (; i + 3 < len; i += 4) {
    uint64x2_t r0 = vld1q_u64(u64 + i);
    uint64x2_t r1 = vld1q_u64(u64 + i + 2);
    r0 = vshrq_n_u64(r0, 12);
    r1 = vshrq_n_u64(r1, 12);
    r0 = vorrq_u64(r0, expo);
    r1 = vorrq_u64(r1, expo);
    float64x2_t d0 = vreinterpretq_f64_u64(r0);
    float64x2_t d1 = vreinterpretq_f64_u64(r1);
    d0 = vfmaq_f64(c, d0, s);
    d1 = vfmaq_f64(c, d1, s);
    vst1q_f64(x + i, d0);
    vst1q_f64(x + i + 2, d1);
  }
#endif
  for (; i < len; i++) {
    uint64_t bits = (u64[i] >> 12) | 0x3ff0000000000000ULL;
    double d;
    memcpy(&d, &bits, sizeof(d));
    x[i] = shift + b*d;
  }
}

ALWAYS_INLINE void rand_dble_4x53(double x[], size_t len, uint64_t u64[]) {
  size_t j = 0;
#if HAVE_NEON
  float64x2_t scale = vdupq_n_f64(0x1.0p-53);
  for (; j + 3 < len; j += 4) {
    uint64x2_t u0 = vld1q_u64(u64 + j);
    uint64x2_t u1 = vld1q_u64(u64 + j + 2);
    u0 = vshrq_n_u64(u0, 11);
    u1 = vshrq_n_u64(u1, 11);
    float64x2_t d0 = vmulq_f64(vcvtq_f64_u64(u0), scale);
    float64x2_t d1 = vmulq_f64(vcvtq_f64_u64(u1), scale);
    vst1q_f64(x + j, d0);
    vst1q_f64(x + j + 2, d1);
  }
#endif
  for (; j < len; j++) {
    x[j] = U64_TO_DOUBLE(u64[j]);
  }
}

ALWAYS_INLINE void rand_dble(double x[], size_t len, randompack_rng *rng) {
#if defined(BUILD_AVX512)
  if (rng->cpu_has_avx512 && !rng->usefullmantissa) {
    rand_dble_avx512(x, len, rng);
    return;
  }
#endif
#if defined(BUILD_AVX2)
  if (rng->cpu_has_avx2) {
    rand_dble_avx2(x, len, rng);
    return;
  }
#endif
  int w = enter_u64_mode(rng);
  size_t i = 0;
  while (i < len) {
    fill_if_empty(rng, &w);
    size_t remain = len - i;
    size_t avail = BUFSIZE - (size_t)w;
    size_t take = remain < avail ? remain : avail;
    if (rng->usefullmantissa) {
      rand_dble_4x53(x + i, take, rng->buf.u64 + (size_t)w);
    }
    else {
      rand_dble_4x52(x + i, take, rng->buf.u64 + (size_t)w);
    }
    w += (int)take;
    i += take;
  }
  exit_u64_mode(rng, w);
}

ALWAYS_INLINE void rand_unif_neon(double x[], size_t len, double a, double b,
  randompack_rng *rng) {
  int w = enter_u64_mode(rng);
  size_t i = 0;
  while (i < len) {
    fill_if_empty(rng, &w);
    size_t remain = len - i;
    size_t avail = BUFSIZE - (size_t)w;
    size_t take = remain < avail ? remain : avail;
    rand_unif_4x52(x + i, take, rng->buf.u64 + (size_t)w, a, b);
    w += (int)take;
    i += take;
  }
  exit_u64_mode(rng, w);
}
