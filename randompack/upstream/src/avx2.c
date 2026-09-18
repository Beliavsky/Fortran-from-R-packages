// -*- C -*-

#include "randompack_internal.h"

#if !defined(__x86_64__) && !defined(_M_X64)

bool cpu_has_avx2(void) {
  return false;
}

void fill_x256ppsimd_avx2(uint64_t *buf, size_t len, randompack_state *state) {
  (void)buf;
  (void)len;
  (void)state;
}

void fill_x256sssimd_avx2(uint64_t *buf, size_t len, randompack_state *state) {
  (void)buf;
  (void)len;
  (void)state;
}

void fill_sfc64simd_avx2(uint64_t *buf, size_t len, randompack_state *state) {
  (void)buf;
  (void)len;
  (void)state;
}

void rand_dble_avx2(double x[], size_t len, randompack_rng *rng) {
  (void)x;
  (void)len;
  (void)rng;
}

void rand_float_avx2(float x[], size_t len, randompack_rng *rng) {
  (void)x;
  (void)len;
  (void)rng;
}

void shift_scale_double_avx2(double x[], size_t len, double shift, double scale) {
  for (size_t i = 0; i < len; i++) x[i] = shift + scale*x[i];
}

void shift_scale_float_avx2(float x[], size_t len, float shift, float scale) {
  for (size_t i = 0; i < len; i++) x[i] = shift + scale*x[i];
}

void scale_double_avx2(double x[], size_t len, double scale) {
  shift_scale_double_avx2(x, len, 0, scale);
}

void scale_float_avx2(float x[], size_t len, float scale) {
  shift_scale_float_avx2(x, len, 0, scale);
}

#else

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <cpuid.h>
#endif

#include <immintrin.h>

#ifndef __m256i_u
typedef __m256i __m256i_u;
#endif

#include <string.h>
#include "buffer_draw.inc"

#if defined(RANDOMPACK_TEST_HOOKS)
#define TEST_HOOK DEFAULT_VISIBLE
static int avx2_used = 0;
TEST_HOOK int randompack_avx2_used(void) {
  return avx2_used;
}
TEST_HOOK void randompack_avx2_reset(void) {
  avx2_used = 0;
}
#endif

HIDDEN bool cpu_has_avx2(void) {
#if defined(_MSC_VER)
  int info[4];
  __cpuid(info, 1);
  int osxsave = (info[2] & (1 << 27)) != 0;
  int avx = (info[2] & (1 << 28)) != 0;
  if (!(osxsave && avx)) return false;
  unsigned long long xcr = _xgetbv(0);
  if ((xcr & 0x6) != 0x6) return false;
  __cpuidex(info, 7, 0);
  return (info[1] & (1 << 5)) != 0;
#elif defined(__GNUC__) || defined(__clang__)
  unsigned int eax, ebx, ecx, edx;
  if (!__get_cpuid(1, &eax, &ebx, &ecx, &edx)) return false;
  if (!(ecx & (1 << 27))) return false;
  if (!(ecx & (1 << 28))) return false;
  unsigned int xcr0_lo, xcr0_hi;
  __asm__ volatile ("xgetbv" : "=a"(xcr0_lo), "=d"(xcr0_hi) : "c"(0));
  if ((xcr0_lo & 0x6) != 0x6) return false;
  if (__get_cpuid_max(0, 0) < 7) return false;
  __cpuid_count(7, 0, eax, ebx, ecx, edx);
  return (ebx & (1 << 5)) != 0;
#else
  return false;
#endif
}

#define VEC_T __m256i
#define VEC_LOAD(p) _mm256_loadu_si256((const __m256i_u *)(p))
#define VEC_STORE(p,x) _mm256_storeu_si256((__m256i_u *)(p),(x))
#define VEC_ADD(a,b) _mm256_add_epi64((a),(b))
#define VEC_XOR(a,b) _mm256_xor_si256((a),(b))
#define VEC_MUL9(a) _mm256_add_epi64(_mm256_slli_epi64((a),3),(a))
#define VEC_SHR11(a) _mm256_srli_epi64((a),11)
#define VEC_SHL17(a) _mm256_slli_epi64((a),17)
#define VEC_ROTL45(x) _mm256_or_si256( \
  _mm256_slli_epi64((x),45), _mm256_srli_epi64((x),64-45))

#define VEC_ROTL23(x) _mm256_or_si256( \
  _mm256_slli_epi64((x),23), _mm256_srli_epi64((x),64-23))
#define VEC_ROTL24(x) _mm256_or_si256( \
  _mm256_slli_epi64((x),24), _mm256_srli_epi64((x),40))

#define FAST_STEP_VEC(s0,s1,s2,s3,outv) do { \
  VEC_T r_; \
  VEC_T t_; \
  r_ = VEC_ADD(VEC_ROTL23(VEC_ADD((s0),(s3))),(s0)); \
  t_ = VEC_SHL17((s1)); \
  (s2) = VEC_XOR((s2),(s0)); \
  (s3) = VEC_XOR((s3),(s1)); \
  (s1) = VEC_XOR((s1),(s2)); \
  (s0) = VEC_XOR((s0),(s3)); \
  (s2) = VEC_XOR((s2),t_); \
  (s3) = VEC_ROTL45((s3)); \
  (outv) = r_; \
} while (0)

// Following 3 macros are used by xoshiro256**
#define VEC_ROTL7(x) _mm256_or_si256( \
  _mm256_slli_epi64((x),7), _mm256_srli_epi64((x),64-7))
#define VEC_MUL5(x) _mm256_add_epi64((x),_mm256_slli_epi64((x),2))
#define FAST_STEP_VEC_SS(s0,s1,s2,s3,outv) do { \
  VEC_T r_; \
  VEC_T t_; \
  r_ = VEC_MUL9(VEC_ROTL7(VEC_MUL5(s1))); \
  t_ = VEC_SHL17((s1)); \
  (s2) = VEC_XOR((s2),(s0)); \
  (s3) = VEC_XOR((s3),(s1)); \
  (s1) = VEC_XOR((s1),(s2)); \
  (s0) = VEC_XOR((s0),(s3)); \
  (s2) = VEC_XOR((s2),t_); \
  (s3) = VEC_ROTL45((s3)); \
  (outv) = r_; \
} while (0)

HIDDEN void fill_x256ppsimd_avx2(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  xo256 *st = &state->xo;
  VEC_T s00 = VEC_LOAD(&st->s0[0]);
  VEC_T s10 = VEC_LOAD(&st->s1[0]);
  VEC_T s20 = VEC_LOAD(&st->s2[0]);
  VEC_T s30 = VEC_LOAD(&st->s3[0]);
  VEC_T s01 = VEC_LOAD(&st->s0[4]);
  VEC_T s11 = VEC_LOAD(&st->s1[4]);
  VEC_T s21 = VEC_LOAD(&st->s2[4]);
  VEC_T s31 = VEC_LOAD(&st->s3[4]);
#if defined(RANDOMPACK_TEST_HOOKS)
  avx2_used++;
#endif
  for (size_t i = 0; i < len; i += 8) {
    VEC_T r0, r1;
    FAST_STEP_VEC(s00, s10, s20, s30, r0);
    FAST_STEP_VEC(s01, s11, s21, s31, r1);
    VEC_STORE(out + i, r0);
    VEC_STORE(out + i + 4, r1);
  }
  VEC_STORE(&st->s0[0], s00);
  VEC_STORE(&st->s1[0], s10);
  VEC_STORE(&st->s2[0], s20);
  VEC_STORE(&st->s3[0], s30);
  VEC_STORE(&st->s0[4], s01);
  VEC_STORE(&st->s1[4], s11);
  VEC_STORE(&st->s2[4], s21);
  VEC_STORE(&st->s3[4], s31);
}

HIDDEN void fill_x256sssimd_avx2(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  xo256 *st = &state->xo;
  VEC_T s00 = VEC_LOAD(&st->s0[0]);
  VEC_T s10 = VEC_LOAD(&st->s1[0]);
  VEC_T s20 = VEC_LOAD(&st->s2[0]);
  VEC_T s30 = VEC_LOAD(&st->s3[0]);
  VEC_T s01 = VEC_LOAD(&st->s0[4]);
  VEC_T s11 = VEC_LOAD(&st->s1[4]);
  VEC_T s21 = VEC_LOAD(&st->s2[4]);
  VEC_T s31 = VEC_LOAD(&st->s3[4]);
  for (size_t i = 0; i < len; i += 8) {
    VEC_T r0, r1;
    FAST_STEP_VEC_SS(s00, s10, s20, s30, r0);
    FAST_STEP_VEC_SS(s01, s11, s21, s31, r1);
    VEC_STORE(out + i, r0);
    VEC_STORE(out + i + 4, r1);
  }
  VEC_STORE(&st->s0[0], s00);
  VEC_STORE(&st->s1[0], s10);
  VEC_STORE(&st->s2[0], s20);
  VEC_STORE(&st->s3[0], s30);
  VEC_STORE(&st->s0[4], s01);
  VEC_STORE(&st->s1[4], s11);
  VEC_STORE(&st->s2[4], s21);
  VEC_STORE(&st->s3[4], s31);
}

#define SFC64_STEP_VEC(a,b,c,ctr,one,outv) do { \
  VEC_T t_; \
  (outv) = VEC_ADD(VEC_ADD((a),(b)),(ctr)); \
  (ctr) = VEC_ADD((ctr),(one)); \
  t_ = VEC_XOR((b), VEC_SHR11((b))); \
  (b) = VEC_MUL9((c)); \
  (c) = VEC_ADD(VEC_ROTL24((c)), (outv)); \
  (a) = t_; \
} while (0)

HIDDEN void fill_sfc64simd_avx2(uint64_t *buf, size_t len, randompack_state *state) {
  uint64_t *out = buf;
  xo256 *st = &state->xo;
  VEC_T a0 = VEC_LOAD(&st->s0[0]);
  VEC_T b0 = VEC_LOAD(&st->s1[0]);
  VEC_T c0 = VEC_LOAD(&st->s2[0]);
  VEC_T d0 = VEC_LOAD(&st->s3[0]);
  VEC_T a1 = VEC_LOAD(&st->s0[4]);
  VEC_T b1 = VEC_LOAD(&st->s1[4]);
  VEC_T c1 = VEC_LOAD(&st->s2[4]);
  VEC_T d1 = VEC_LOAD(&st->s3[4]);
  VEC_T one = _mm256_set1_epi64x(1);
#if defined(RANDOMPACK_TEST_HOOKS)
  avx2_used++;
#endif
  for (size_t i = 0; i < len; i += 8) {
    VEC_T r0, r1;
    SFC64_STEP_VEC(a0, b0, c0, d0, one, r0);
    SFC64_STEP_VEC(a1, b1, c1, d1, one, r1);
    VEC_STORE(out + i, r0);
    VEC_STORE(out + i + 4, r1);
  }
  VEC_STORE(&st->s0[0], a0);
  VEC_STORE(&st->s1[0], b0);
  VEC_STORE(&st->s2[0], c0);
  VEC_STORE(&st->s3[0], d0);
  VEC_STORE(&st->s0[4], a1);
  VEC_STORE(&st->s1[4], b1);
  VEC_STORE(&st->s2[4], c1);
  VEC_STORE(&st->s3[4], d1);
}

HIDDEN void rand_dble_avx2(double x[], size_t len, randompack_rng *rng) {
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
      __m256i expo = _mm256_set1_epi64x(0x3ff0000000000000LL);
      __m256d one = _mm256_set1_pd(1);
      for (; j + 3 < take; j += 4) {
        __m256i r = _mm256_loadu_si256((const __m256i_u *)(u64 + j));
        r = _mm256_srli_epi64(r, 12);
        r = _mm256_or_si256(r, expo);
        __m256d d = _mm256_castsi256_pd(r);
        d = _mm256_sub_pd(d, one);
        _mm256_storeu_pd(x + i + j, d);
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

HIDDEN void rand_unif_avx2(double x[], size_t len, double a, double b,
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
    __m256i expo = _mm256_set1_epi64x(0x3ff0000000000000LL);
    __m256d s = _mm256_set1_pd(scale);
    __m256d c = _mm256_set1_pd(shift);
    for (; j + 3 < take; j += 4) {
      __m256i r = _mm256_loadu_si256((const __m256i_u *)(u64 + j));
      r = _mm256_srli_epi64(r, 12);
      r = _mm256_or_si256(r, expo);
      __m256d d = _mm256_castsi256_pd(r);
#if defined(__FMA__)
      d = _mm256_fmadd_pd(d, s, c);
#else
      d = _mm256_add_pd(_mm256_mul_pd(d, s), c);
#endif
      _mm256_storeu_pd(x + i + j, d);
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

HIDDEN void rand_float_avx2(float x[], size_t len, randompack_rng *rng) {
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
      __m256 scale = _mm256_set1_ps(0x1.0p-24f);
      for (; j + 7 < take; j += 8) {
        __m256i r = _mm256_loadu_si256((const __m256i_u *)(u32 + j));
        r = _mm256_srli_epi32(r, 8);
        __m256 d = _mm256_cvtepi32_ps(r);
        d = _mm256_mul_ps(d, scale);
        _mm256_storeu_ps(x + i + j, d);
      }
      for (; j < take; j++) {
        x[i + j] = U32_TO_FLOAT(u32[j]);
      }
    }
    else {
      size_t j = 0;
      __m256 scale = _mm256_set1_ps(0x1.0p-23f);
      for (; j + 7 < take; j += 8) {
        __m256i r = _mm256_loadu_si256((const __m256i_u *)(u32 + j));
        r = _mm256_srli_epi32(r, 9);
        __m256 d = _mm256_cvtepi32_ps(r);
        d = _mm256_mul_ps(d, scale);
        _mm256_storeu_ps(x + i + j, d);
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

HIDDEN void shift_scale_double_avx2(double x[], size_t len, double shift,
  double scale) {
#if defined(RANDOMPACK_TEST_HOOKS)
  avx2_used++;
#endif
  size_t i = 0;
  __m256d s = _mm256_set1_pd(scale);
  __m256d b = _mm256_set1_pd(shift);
  for (; i + 4 <= len; i += 4) {
    __m256d v = _mm256_loadu_pd(x + i);
#if defined(__FMA__)
    v = _mm256_fmadd_pd(v, s, b);
#else
    v = _mm256_add_pd(_mm256_mul_pd(v, s), b);
#endif
    _mm256_storeu_pd(x + i, v);
  }
  for (; i < len; i++) x[i] = shift + scale*x[i];
}

HIDDEN void shift_scale_float_avx2(float x[], size_t len, float shift,
  float scale) {
#if defined(RANDOMPACK_TEST_HOOKS)
  avx2_used++;
#endif
  size_t i = 0;
  __m256 s = _mm256_set1_ps(scale);
  __m256 b = _mm256_set1_ps(shift);
  for (; i + 8 <= len; i += 8) {
    __m256 v = _mm256_loadu_ps(x + i);
#if defined(__FMA__)
    v = _mm256_fmadd_ps(v, s, b);
#else
    v = _mm256_add_ps(_mm256_mul_ps(v, s), b);
#endif
    _mm256_storeu_ps(x + i, v);
  }
  for (; i < len; i++) x[i] = shift + scale*x[i];
}

HIDDEN void scale_double_avx2(double x[], size_t len, double scale) {
  shift_scale_double_avx2(x, len, 0, scale);
}

HIDDEN void scale_float_avx2(float x[], size_t len, float scale) {
  shift_scale_float_avx2(x, len, 0, scale);
}

#endif
