// -*- C -*-

#if HAVE_NEON
#include <arm_neon.h>
#endif

ALWAYS_INLINE void shift_scale_double_inplace(double *x, size_t len, double shift,
  double scale, randompack_rng *rng) {
  if (!rng->bitexact) {
#if defined(BUILD_AVX512)
    if (rng->cpu_has_avx512) {
      shift_scale_double_avx512(x, len, shift, scale);
      return;
    }
#endif
#if defined(BUILD_AVX2)
    if (rng->cpu_has_avx2) {
      shift_scale_double_avx2(x, len, shift, scale);
      return;
    }
#endif
#if HAVE_NEON && defined(__aarch64__)
    size_t i = 0;
    float64x2_t s = vdupq_n_f64(scale);
    float64x2_t b = vdupq_n_f64(shift);
    for (; i + 4 <= len; i += 4) {
      float64x2_t v0 = vld1q_f64(x + i);
      float64x2_t v1 = vld1q_f64(x + i + 2);
      v0 = vfmaq_f64(b, v0, s);
      v1 = vfmaq_f64(b, v1, s);
      vst1q_f64(x + i, v0);
      vst1q_f64(x + i + 2, v1);
    }
    for (; i + 2 <= len; i += 2) {
      float64x2_t v = vld1q_f64(x + i);
      v = vfmaq_f64(b, v, s);
      vst1q_f64(x + i, v);
    }
    for (; i < len; i++) x[i] = fma(scale, x[i], shift);
    return;
#endif
  }
  (void)rng;
  for (size_t i = 0; i < len; i++) x[i] = shift + scale*x[i];
}

static void scale_double_inplace(double *x, size_t len, double scale,
  randompack_rng *rng) {
  shift_scale_double_inplace(x, len, 0, scale, rng);
}
