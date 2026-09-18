// -*- C -*-
// Sleef in-place log/exp helpers.

#include <math.h>
#include <stddef.h>
#include "randompack_config.h"

#if defined(__x86_64__) || defined(_M_X64)

#include <immintrin.h>

void sleef_expd4_u10avx2_array(double *x, size_t len);
void sleef_logd4_u35avx2_array(double *x, size_t len);
void sleef_expf8_u10avx2_array(float *x, size_t len);
void sleef_logf8_u10avx2_array(float *x, size_t len);

HIDDEN void sleef_exp_inplace(double *x, size_t len) {
  sleef_expd4_u10avx2_array(x, len);
}

HIDDEN void sleef_log_inplace(double *x, size_t len) {
  sleef_logd4_u35avx2_array(x, len);
}

HIDDEN void sleef_expf_inplace(float *x, size_t len) {
  sleef_expf8_u10avx2_array(x, len);
}

HIDDEN void sleef_logf_inplace(float *x, size_t len) {
  sleef_logf8_u10avx2_array(x, len);
}

#elif (defined(__aarch64__) || defined(_M_ARM64)) && !defined(__APPLE__)

#include <arm_neon.h>

float64x2_t Sleef_expd2_u10advsimd(float64x2_t d);
float64x2_t Sleef_logd2_u35advsimd(float64x2_t d);
float32x4_t Sleef_expf4_u10advsimd(float32x4_t d);
float32x4_t Sleef_logf4_u10advsimd(float32x4_t d);

HIDDEN void sleef_exp_inplace(double *x, size_t len) {
  size_t i = 0;
  for (; i + 2 <= len; i += 2) {
    float64x2_t d = vld1q_f64(x + i);
    d = Sleef_expd2_u10advsimd(d);
    vst1q_f64(x + i, d);
  }
  for (; i < len; i++) x[i] = exp(x[i]);
}

HIDDEN void sleef_log_inplace(double *x, size_t len) {
  size_t i = 0;
  for (; i + 2 <= len; i += 2) {
    float64x2_t d = vld1q_f64(x + i);
    d = Sleef_logd2_u35advsimd(d);
    vst1q_f64(x + i, d);
  }
  for (; i < len; i++) x[i] = log(x[i]);
}

HIDDEN void sleef_expf_inplace(float *x, size_t len) {
  size_t i = 0;
  for (; i + 4 <= len; i += 4) {
    float32x4_t d = vld1q_f32(x + i);
    d = Sleef_expf4_u10advsimd(d);
    vst1q_f32(x + i, d);
  }
  for (; i < len; i++) x[i] = expf(x[i]);
}

HIDDEN void sleef_logf_inplace(float *x, size_t len) {
  size_t i = 0;
  for (; i + 4 <= len; i += 4) {
    float32x4_t d = vld1q_f32(x + i);
    d = Sleef_logf4_u10advsimd(d);
    vst1q_f32(x + i, d);
  }
  for (; i < len; i++) x[i] = logf(x[i]);
}

#else

HIDDEN void sleef_exp_inplace(double *x, size_t len) {
  for (size_t i = 0; i < len; i++) x[i] = exp(x[i]);
}

HIDDEN void sleef_log_inplace(double *x, size_t len) {
  for (size_t i = 0; i < len; i++) x[i] = log(x[i]);
}

HIDDEN void sleef_expf_inplace(float *x, size_t len) {
  for (size_t i = 0; i < len; i++) x[i] = expf(x[i]);
}

HIDDEN void sleef_logf_inplace(float *x, size_t len) {
  for (size_t i = 0; i < len; i++) x[i] = logf(x[i]);
}
#endif
