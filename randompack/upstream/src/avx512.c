// -*- C -*-

#include "randompack_internal.h"

#if !defined(__x86_64__) && !defined(_M_X64)

bool cpu_has_avx512(void) {
  return false;
}

void fill_x256ppsimd_avx512(uint64_t *buf, size_t len,
  randompack_state *state) {
  (void)buf;
  (void)len;
  (void)state;
}

void fill_x256sssimd_avx512(uint64_t *buf, size_t len, randompack_state *state) {
  (void)buf;
  (void)len;
  (void)state;
}

void fill_sfc64simd_avx512(uint64_t *buf, size_t len, randompack_state *state) {
  (void)buf;
  (void)len;
  (void)state;
}

void rand_dble_avx512(double x[], size_t len, randompack_rng *rng) {
  (void)x;
  (void)len;
  (void)rng;
}

void rand_float_avx512(float x[], size_t len, randompack_rng *rng) {
  (void)x;
  (void)len;
  (void)rng;
}

void shift_scale_double_avx512(double x[], size_t len, double shift, double scale) {
  for (size_t i = 0; i < len; i++) x[i] = shift + scale*x[i];
}

void shift_scale_float_avx512(float x[], size_t len, float shift, float scale) {
  for (size_t i = 0; i < len; i++) x[i] = shift + scale*x[i];
}

void scale_double_avx512(double x[], size_t len, double scale) {
  shift_scale_double_avx512(x, len, 0, scale);
}

void scale_float_avx512(float x[], size_t len, float scale) {
  shift_scale_float_avx512(x, len, 0, scale);
}

#else

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <cpuid.h>
#endif

#include <immintrin.h>

#ifndef __m512i_u
typedef __m512i __m512i_u;
#endif

#include <string.h>
#include "buffer_draw.inc"

HIDDEN bool cpu_has_avx512(void) {
#if defined(_MSC_VER)
  int info[4];
  __cpuid(info, 1);
  int osxsave = (info[2] & (1 << 27)) != 0;
  int avx = (info[2] & (1 << 28)) != 0;
  if (!(osxsave && avx)) return false;
  unsigned long long xcr = _xgetbv(0);
  if ((xcr & 0xe6) != 0xe6) return false;
  __cpuidex(info, 7, 0);
  return (info[1] & (1 << 16)) != 0 && (info[1] & (1 << 17)) != 0;
#elif defined(__GNUC__) || defined(__clang__)
  unsigned int eax, ebx, ecx, edx;
  if (!__get_cpuid(1, &eax, &ebx, &ecx, &edx)) return false;
  if (!(ecx & (1 << 27))) return false;
  if (!(ecx & (1 << 28))) return false;
  unsigned int xcr0_lo, xcr0_hi;
  __asm__ volatile ("xgetbv" : "=a"(xcr0_lo), "=d"(xcr0_hi) : "c"(0));
  if ((xcr0_lo & 0xe6) != 0xe6) return false;
  if (__get_cpuid_max(0, 0) < 7) return false;
  __cpuid_count(7, 0, eax, ebx, ecx, edx);
  return (ebx & (1 << 16)) != 0 && (ebx & (1 << 17)) != 0;
#else
  return false;
#endif
}

HIDDEN void fill_x256ppsimd_avx512(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  xo256 *st = &state->xo;
  __m512i s0 = _mm512_loadu_si512((const void *)&st->s0[0]);
  __m512i s1 = _mm512_loadu_si512((const void *)&st->s1[0]);
  __m512i s2 = _mm512_loadu_si512((const void *)&st->s2[0]);
  __m512i s3 = _mm512_loadu_si512((const void *)&st->s3[0]);
  size_t i = 0;
  for (; i < len; i += 8) {
    __m512i sum0 = _mm512_add_epi64(s0, s3);
    __m512i r0 = _mm512_add_epi64(
      _mm512_or_si512(_mm512_slli_epi64(sum0, 23), _mm512_srli_epi64(sum0, 41)),
      s0
    );
    __m512i t0 = _mm512_slli_epi64(s1, 17);
    s2 = _mm512_xor_si512(s2, s0);
    s3 = _mm512_xor_si512(s3, s1);
    s1 = _mm512_xor_si512(s1, s2);
    s0 = _mm512_xor_si512(s0, s3);
    s2 = _mm512_xor_si512(s2, t0);
    s3 = _mm512_or_si512(_mm512_slli_epi64(s3, 45), _mm512_srli_epi64(s3, 19));
    _mm512_storeu_si512((void *)(out + i), r0);
  }
  _mm512_storeu_si512((void *)&st->s0[0], s0);
  _mm512_storeu_si512((void *)&st->s1[0], s1);
  _mm512_storeu_si512((void *)&st->s2[0], s2);
  _mm512_storeu_si512((void *)&st->s3[0], s3);
}

HIDDEN void fill_x256sssimd_avx512(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  xo256 *st = &state->xo;
  __m512i s0 = _mm512_loadu_si512((const void *)&st->s0[0]);
  __m512i s1 = _mm512_loadu_si512((const void *)&st->s1[0]);
  __m512i s2 = _mm512_loadu_si512((const void *)&st->s2[0]);
  __m512i s3 = _mm512_loadu_si512((const void *)&st->s3[0]);
  size_t i = 0;
  for (; i < len; i += 8) {
    __m512i r0 = _mm512_add_epi64(_mm512_slli_epi64(s1, 2), s1);
    r0 = _mm512_or_si512(_mm512_slli_epi64(r0, 7), _mm512_srli_epi64(r0, 57));
    r0 = _mm512_add_epi64(_mm512_slli_epi64(r0, 3), r0);
    __m512i t0 = _mm512_slli_epi64(s1, 17);
    s2 = _mm512_xor_si512(s2, s0);
    s3 = _mm512_xor_si512(s3, s1);
    s1 = _mm512_xor_si512(s1, s2);
    s0 = _mm512_xor_si512(s0, s3);
    s2 = _mm512_xor_si512(s2, t0);
    s3 = _mm512_or_si512(_mm512_slli_epi64(s3, 45), _mm512_srli_epi64(s3, 19));
    _mm512_storeu_si512((void *)(out + i), r0);
  }
  _mm512_storeu_si512((void *)&st->s0[0], s0);
  _mm512_storeu_si512((void *)&st->s1[0], s1);
  _mm512_storeu_si512((void *)&st->s2[0], s2);
  _mm512_storeu_si512((void *)&st->s3[0], s3);
}

HIDDEN void fill_sfc64simd_avx512(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  xo256 *st = &state->xo;
  __m512i a = _mm512_loadu_si512((const void *)&st->s0[0]);
  __m512i b = _mm512_loadu_si512((const void *)&st->s1[0]);
  __m512i c = _mm512_loadu_si512((const void *)&st->s2[0]);
  __m512i d = _mm512_loadu_si512((const void *)&st->s3[0]);
  __m512i one = _mm512_set1_epi64(1);
  size_t i = 0;
  for (; i < len; i += 8) {
    __m512i t0;
    __m512i r0 = _mm512_add_epi64(_mm512_add_epi64(a, b), d);
    d = _mm512_add_epi64(d, one);
    t0 = _mm512_xor_si512(b, _mm512_srli_epi64(b, 11));
    b = _mm512_add_epi64(c, _mm512_slli_epi64(c, 3));
    c = _mm512_add_epi64(
      _mm512_or_si512(_mm512_slli_epi64(c, 24), _mm512_srli_epi64(c, 40)),
      r0
    );
    a = t0;
    _mm512_storeu_si512((void *)(out + i), r0);
  }
  _mm512_storeu_si512((void *)&st->s0[0], a);
  _mm512_storeu_si512((void *)&st->s1[0], b);
  _mm512_storeu_si512((void *)&st->s2[0], c);
  _mm512_storeu_si512((void *)&st->s3[0], d);
}

HIDDEN void rand_float_avx512(float x[], size_t len, randompack_rng *rng) {
  int w = enter_u32_mode(rng);
  size_t i = 0;
  while (i < len) {
    fill_if_empty32(rng, &w);
    size_t remain = len - i;
    size_t avail = 2*BUFSIZE - (size_t)w;
    size_t take = remain < avail ? remain : avail;
    uint32_t *u32 = rng->buf.u32 + (size_t)w;
    if (rng->usefullmantissa) {
      size_t j = 0;
      __m512 scale = _mm512_set1_ps(0x1.0p-24f);
      for (; j + 15 < take; j += 16) {
        __m512i r = _mm512_loadu_si512((const void *)(u32 + j));
        r = _mm512_srli_epi32(r, 8);
        __m512 d = _mm512_cvtepi32_ps(r);
        d = _mm512_mul_ps(d, scale);
        _mm512_storeu_ps(x + i + j, d);
      }
      for (; j < take; j++) {
        x[i + j] = U32_TO_FLOAT(u32[j]);
      }
    }
    else {
      size_t j = 0;
      __m512 scale = _mm512_set1_ps(0x1.0p-23f);
      for (; j + 15 < take; j += 16) {
        __m512i r = _mm512_loadu_si512((const void *)(u32 + j));
        r = _mm512_srli_epi32(r, 9);
        __m512 d = _mm512_cvtepi32_ps(r);
        d = _mm512_mul_ps(d, scale);
        _mm512_storeu_ps(x + i + j, d);
      }
      for (; j < take; j++) {
        x[i + j] = (float)(u32[j] >> 9)*0x1.0p-23f;
      }
    }
    w += (int)take;
    i += take;
  }
  exit_u32_mode(rng, w);
}

HIDDEN void shift_scale_double_avx512(double x[], size_t len, double shift,
  double scale) {
  size_t i = 0;
  __m512d s = _mm512_set1_pd(scale);
  __m512d b = _mm512_set1_pd(shift);
  for (; i + 8 <= len; i += 8) {
    __m512d v = _mm512_loadu_pd(x + i);
#if defined(__FMA__)
    v = _mm512_fmadd_pd(v, s, b);
#else
    v = _mm512_add_pd(_mm512_mul_pd(v, s), b);
#endif
    _mm512_storeu_pd(x + i, v);
  }
  for (; i < len; i++) x[i] = shift + scale*x[i];
}

HIDDEN void shift_scale_float_avx512(float x[], size_t len, float shift,
  float scale) {
  size_t i = 0;
  __m512 s = _mm512_set1_ps(scale);
  __m512 b = _mm512_set1_ps(shift);
  for (; i + 16 <= len; i += 16) {
    __m512 v = _mm512_loadu_ps(x + i);
#if defined(__FMA__)
    v = _mm512_fmadd_ps(v, s, b);
#else
    v = _mm512_add_ps(_mm512_mul_ps(v, s), b);
#endif
    _mm512_storeu_ps(x + i, v);
  }
  for (; i < len; i++) x[i] = shift + scale*x[i];
}

HIDDEN void scale_double_avx512(double x[], size_t len, double scale) {
  shift_scale_double_avx512(x, len, 0, scale);
}

HIDDEN void scale_float_avx512(float x[], size_t len, float scale) {
  shift_scale_float_avx512(x, len, 0, scale);
}

HIDDEN void rand_dble_avx512(double x[], size_t len, randompack_rng *rng) {
  int w = enter_u64_mode(rng);
  size_t i = 0;
  while (i < len) {
    fill_if_empty(rng, &w);
    size_t remain = len - i;
    size_t avail = BUFSIZE - (size_t)w;
    size_t take = remain < avail ? remain : avail;
    uint64_t *u64 = rng->buf.u64 + (size_t)w;
    if (rng->usefullmantissa) {
      for (size_t j = 0; j < take; j++) x[i + j] = U64_TO_DOUBLE(u64[j]);
    }
    else {
      size_t j = 0;
      __m512i expo = _mm512_set1_epi64(0x3ff0000000000000ULL);
      __m512d one = _mm512_set1_pd(1);
      for (; j + 7 < take; j += 8) {
        __m512i r = _mm512_loadu_si512((const void *)(u64 + j));
        r = _mm512_srli_epi64(r, 12);
        r = _mm512_or_si512(r, expo);
        __m512d d = _mm512_castsi512_pd(r);
        d = _mm512_sub_pd(d, one);
        _mm512_storeu_pd(x + i + j, d);
      }
      for (; j < take; j++) {
        uint64_t bits = (u64[j] >> 12) | 0x3ff0000000000000ULL;
        double d;
        memcpy(&d, &bits, sizeof(d));
        x[i + j] = d - 1;
      }
    }
    w += (int)take;
    i += take;
  }
  exit_u64_mode(rng, w);
}

HIDDEN void rand_unif_avx512(double x[], size_t len, double a, double b,
  randompack_rng *rng) {
  int w = enter_u64_mode(rng);
  double scale = b;
  double shift = a - b;
  size_t i = 0;
  while (i < len) {
    fill_if_empty(rng, &w);
    size_t remain = len - i;
    size_t avail = BUFSIZE - (size_t)w;
    size_t take = remain < avail ? remain : avail;
    uint64_t *u64 = rng->buf.u64 + (size_t)w;
    size_t j = 0;
    __m512i expo = _mm512_set1_epi64(0x3ff0000000000000ULL);
    __m512d s = _mm512_set1_pd(scale);
    __m512d c = _mm512_set1_pd(shift);
    for (; j + 7 < take; j += 8) {
      __m512i r = _mm512_loadu_si512((const void *)(u64 + j));
      r = _mm512_srli_epi64(r, 12);
      r = _mm512_or_si512(r, expo);
      __m512d d = _mm512_castsi512_pd(r);
#if defined(__FMA__)
      d = _mm512_fmadd_pd(d, s, c);
#else
      d = _mm512_add_pd(_mm512_mul_pd(d, s), c);
#endif
      _mm512_storeu_pd(x + i + j, d);
    }
    for (; j < take; j++) {
      uint64_t bits = (u64[j] >> 12) | 0x3ff0000000000000ULL;
      double d;
      memcpy(&d, &bits, sizeof(d));
      x[i + j] = shift + scale*d;
    }
    w += (int)take;
    i += take;
  }
  exit_u64_mode(rng, w);
}

#endif
