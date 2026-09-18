// -*- c -*-
#if defined(USE_ACCEL_VV)
void vvexp(double *y, const double *x, const int *n);
void vvlog(double *y, const double *x, const int *n);
#endif

#if defined(BUILD_AVX2) || (defined(__aarch64__) && !defined(__APPLE__))
void sleef_exp_inplace(double *x, size_t len);
void sleef_log_inplace(double *x, size_t len);
#endif
#if defined(BUILD_AVX512)
void sleef_exp_inplace_avx512(double *x, size_t len);
void sleef_log_inplace_avx512(double *x, size_t len);
#endif

NEVER_INLINE void exp_inplace(double *x, size_t len, randompack_rng *rng) {
  if (!rng->bitexact) {
#if defined(USE_ACCEL_VV)
    for (size_t i = 0; i < len; ) {
      size_t chunk = len - i > (size_t)INT_MAX ? (size_t)INT_MAX : len - i;
      int n = (int)chunk;
      vvexp(x + i, x + i, &n);
      i += chunk;
    }
    return;
#elif defined(__x86_64__) || defined(_M_X64)
#if defined(BUILD_AVX512)
    if (rng->cpu_has_avx512) {
      sleef_exp_inplace_avx512(x, len);
      return;
    }
#endif
#if defined(BUILD_AVX2)
    if (rng->cpu_has_avx2) {
      sleef_exp_inplace(x, len);
      return;
    }
#endif
    for (size_t i = 0; i < len; i++) x[i] = exp(x[i]);
    return;
// #elif defined(__aarch64__) && !defined(__APPLE__)
//     sleef_exp_inplace(x, len);
//     return;
#else
    for (size_t i = 0; i < len; i++) x[i] = exp(x[i]);
    return;
#endif
  }
  for (size_t i = 0; i < len; i++) x[i] = openlibm_exp(x[i]);
}

NEVER_INLINE void log_inplace(double *x, size_t len, randompack_rng *rng) {
  if (!rng->bitexact) {
#if defined(USE_ACCEL_VV)
    for (size_t i = 0; i < len; ) {
      size_t chunk = len - i > (size_t)INT_MAX ? (size_t)INT_MAX : len - i;
      int n = (int)chunk;
      vvlog(x + i, x + i, &n);
      i += chunk;
    }
    return;
#elif defined(__x86_64__) || defined(_M_X64)
#if defined(BUILD_AVX512)
    if (rng->cpu_has_avx512) {
      sleef_log_inplace_avx512(x, len);
      return;
    }
#endif
#if defined(BUILD_AVX2)
    if (rng->cpu_has_avx2) {
      sleef_log_inplace(x, len);
      return;
    }
#endif
    for (size_t i = 0; i < len; i++) x[i] = log(x[i]);
    return;
// #elif defined(__aarch64__) && !defined(__APPLE__)
//     sleef_log_inplace(x, len);
//     return;
#else
    for (size_t i = 0; i < len; i++) x[i] = log(x[i]);
    return;
#endif
  }
  for (size_t i = 0; i < len; i++) x[i] = openlibm_log(x[i]);
}
