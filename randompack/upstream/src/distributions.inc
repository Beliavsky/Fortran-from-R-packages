// -*- mode: c; -*-
#include <float.h>
#if defined(USE_ACCEL_VV)
void vvexp(double *y, const double *x, const int *n);
void vvlog(double *y, const double *x, const int *n);
#endif

static int LastNonzeroColumn(int m, int n, double A[]) {
  // If they are all zero, return 0.
  int k = n-1;
  while (k >= 0) {
    int imx = iamax(m, A+k*m, 1);
    if (A[imx + k*m] > 0) return k + 1;
    k--;
  }
  return k + 2;
}

//#include "zmgd.inc"   // static zmgd() lives here
#include "norm_exp.inc"
#include "rand_dble.inc"

static void rand_lognormal(double x[], size_t len, double mu, double sigma,
  randompack_rng *rng) {
  rand_norm(x, len, rng); // N(0,1)
  if (mu != 0 || sigma != 1) shift_scale_double_inplace(x, len, mu, sigma, rng);
  exp_inplace(x, len, rng);
}

// All the rand_uintXX, rand_int and rand_long_long are strict-aliasing-safe
static void rand_uint8(uint8_t x[], size_t len, uint8_t bound, randompack_rng *rng) {
  size_t i = 0;
  uint32_t thresh = 0x100u % bound;
  while (i < len) {
    uint32_t r = draw_byte(rng);
    uint32_t prod = r*bound;
    uint32_t lo = prod & 0xffu;
    if (lo >= thresh) x[i++] = (uint8_t)(prod >> 8);
  }
}

static void rand_uint16(uint16_t x[], size_t len, uint16_t bound, randompack_rng *rng) {
  size_t i = 0;
  uint32_t thresh = (uint32_t)((uint32_t)0x10000u % bound);
  int w = enter_u16_mode(rng);
  while (i < len) {
    uint32_t r = draw_u16(rng, &w);
    uint32_t prod = r*bound;
    uint32_t lo = prod & 0xffffu;
    if (lo >= thresh) x[i++] = (uint16_t)(prod >> 16);
  }
  exit_u16_mode(rng, w);
}

static void rand_uint32(uint32_t x[], size_t len, uint32_t bound, randompack_rng *rng) {
  size_t i = 0;
  if (bound <= 1000u) {
    uint32_t thresh = 0x10000u % bound;
    int w = enter_u16_mode(rng);
    while (i < len) {
      uint32_t r = draw_u16(rng, &w);
      uint32_t prod = r*bound;
      uint32_t lo = prod & 0xffffu;
      if (lo >= thresh) x[i++] = prod >> 16;
    }
    exit_u16_mode(rng, w);
  }
  else {
    uint32_t thresh = ((~bound) + 1u) % bound;
    int w = enter_u32_mode(rng);
    while (i < len) {
      uint32_t r = draw_u32(rng, &w);
      uint64_t prod = (uint64_t)r*(uint64_t)bound;
      uint32_t lo = (uint32_t)prod;
      if (lo >= thresh) x[i++] = (uint32_t)(prod >> 32);
    }
    exit_u32_mode(rng, w);
  }
}

static void rand_uint64(uint64_t x[], size_t len, uint64_t bound, randompack_rng *rng) {
  size_t i = 0;
  if (bound <= 10000ull) {
    uint64_t thresh = 0x10000ull % bound;
    int w = enter_u16_mode(rng);
    while (i < len) {
      uint64_t r = draw_u16(rng, &w);
      uint64_t prod = r*bound;
      uint64_t lo = prod & 0xffffull;
      if (lo >= thresh) x[i++] = prod >> 16;
    }
    exit_u16_mode(rng, w);
  }
  else {
    uint64_t thresh = ((~bound) + 1ull) % bound;
    int w = enter_u64_mode(rng);
    while (i < len) {
      uint64_t r = draw_u64(rng, &w);
      uint64_t hi, lo;
      MUL64_WIDE(r, bound, hi, lo);
      if (lo >= thresh) x[i++] = hi;
    }
    exit_u64_mode(rng, w);
  }
}

static void rand_int(int x[], size_t len, int m, uint32_t bound, randompack_rng *rng) {
  size_t i = 0;
  if (bound <= 1000) {
    uint32_t thresh = 0x10000u % bound;
    int w = enter_u16_mode(rng);
    while (i < len) {
      uint32_t r = draw_u16(rng, &w);
      uint32_t prod = r*bound;
      uint32_t lo = prod & 0xffffu;
      if (lo >= thresh) x[i++] = m + (int)(prod >> 16);
    }
    exit_u16_mode(rng, w);
  }
  else {
    uint32_t thresh = ((~bound) + 1u) % bound;
    int w = enter_u32_mode(rng);
    while (i < len) {
      uint32_t r = draw_u32(rng, &w);
      uint64_t prod = (uint64_t)r*(uint64_t)bound;
      uint32_t lo = (uint32_t)prod;
      if (lo >= thresh) x[i++] = m + (int)(prod >> 32);
    }
    exit_u32_mode(rng, w);
  }
}

static void rand_long_long(long long x[], size_t len, long long m, uint64_t bound,
                           randompack_rng *rng) {
  size_t i = 0;
  if (bound <= 10000ull) {
    uint64_t thresh = 0x10000ull % bound;
    int w = enter_u16_mode(rng);
    while (i < len) {
      uint64_t r = draw_u16(rng, &w);
      uint64_t prod = r*bound;
      uint64_t lo = prod & 0xffffull;
      if (lo >= thresh) x[i++] = m + (long long)(prod >> 16);
    }
    exit_u16_mode(rng, w);
  }
  else {
    uint64_t thresh = ((~bound) + 1ull) % bound;
    int w = enter_u64_mode(rng);
    while (i < len) {
      uint64_t r = draw_u64(rng, &w);
      uint64_t hi, lo;
      MUL64_WIDE(r, bound, hi, lo);
      if (lo >= thresh) x[i++] = m + (long long)hi;
    }
    exit_u64_mode(rng, w);
  }
}

static void rand_perm(int p[], int len, randompack_rng *rng) {
  // Fisher-Yates algorithm, stores a random permutation of 0..len-1 in p
  for (int i = 0; i < len; i++) {
    int j;
    rand_int(&j, 1, 0, i + 1, rng);
    if (j != i) p[i] = p[j];
    p[j] = i;
  }
}

static inline bool hashtable_has(int *keys, int mask, int v) {
  uint32_t h = mix32((uint32_t)v);
  int i = (int)(h & (uint32_t)mask);
  for (;;) {
    int k = keys[i];
    if (k == -1) return false;
    if (k == v) return true;
    i = (i + 1) & mask;
  }
}

static inline void hashtable_add(int *keys, int mask, int v) {
  uint32_t h = mix32((uint32_t)v);
  int i = (int)(h & (uint32_t)mask);
  for (;;) {
    int k = keys[i];
    if (k == -1 || k == v) { keys[i] = v; return; }
    i = (i + 1) & mask;
  }
}

static void rand_sample_reservoir(int p[], int len, int k, randompack_rng *rng) {
  // Select k values from 0..len-1, c.f. reservoir sampling on Wikipedia
  int r;
  p[0] = 0;
  for (int i = 0; i < k; i++) {
    rand_int(&r, 1, 0, i + 1, rng);
    if (r != i) p[i] = p[r];
    p[r] = i;
  }
  for (int i = k; i < len; i++) {
    rand_int(&r, 1, 0, i + 1, rng);
    if (r < k) p[r] = i;
  }
}

static void rand_sample_floyd(int p[], int len, int k, randompack_rng *rng) {
  // Floyd's algorithm: O(k) expected. p[0..k-1] gets k distinct integers from 0..len-1.
  if (k <= 0) return;
  int m = 1;
  while (m < 2*k) m <<= 1;
  int mask = m - 1;
  int *keys;
  if (!ALLOC(keys, m)) { // fall-back
    rand_sample_reservoir(p, len, k, rng);
    return;
  }
  for (int i = 0; i < m; i++) keys[i] = -1;
  int out = 0;
  for (int j = len - k; j < len; j++) {
    int t;
    rand_int(&t, 1, 0, j + 1, rng);
    int pick = hashtable_has(keys, mask, t) ? j : t;
    hashtable_add(keys, mask, pick);
    p[out++] = pick;
  }
  FREE(keys);
}

static void rand_sample(int p[], int len, int k, randompack_rng *rng) {
  // Floyd for small k, reservoir for big k
  if (k < len/2) rand_sample_floyd(p, len, k, rng);
  else rand_sample_reservoir(p, len, k, rng);
}

// -*- C -*-
// Gamma RNG (shape, scale): Marsaglia–Tsang for shape>=1, boost trick for shape<1.

#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include "randompack.h"
#include "randompack_config.h"

ALWAYS_INLINE void gen_gamma_shape1(double *x, size_t len, double shape,
												double scale, randompack_rng *rng) {
  double d = shape - 1/3.0, dscale = d*scale;
  double c = 1/sqrt(9*d), v3, z, v;
  rand_norm(x, len, rng);
  int w = enter_u64_mode(rng);
  for (size_t i = 0; i < len; i++) {
    for (;;) {
      z = x[i];
      v = 1 + c*z;
      if (v > 0) {
        double u, z2 = z*z, z4 = z2*z2, logu, logv3;
        v3 = v*v*v;
        u = draw_u01(rng, &w);
        if (u < 1 - 0.0331*z4) goto accept;
        // Use fixed logs here because this test controls RNG consumption.
        logu = openlibm_log(u);
        logv3 = openlibm_log(v3);
        if (logu < 0.5*z2 + d*(1 - v3 + logv3)) goto accept;
      }
      x[i] = draw_norm(rng, &w);
      continue;
    accept:
      x[i] = dscale*v3;
      break;
    }    
  }
  exit_u64_mode(rng, w);
}

ALWAYS_INLINE void gen_gamma(double *x, size_t len, double shape, double scale,
  randompack_rng *rng) {
  if (shape >= 1) {
    gen_gamma_shape1(x, len, shape, scale, rng);
  }
  else {  // shape < 1
    double inv_shape = 1/shape;
    double t[128];
    gen_gamma_shape1(x, len, shape + 1, scale, rng);
    size_t i = 0;
    while (i < len) {
      size_t n = len - i < 128 ? len - i : 128;
      rand_dble(t, n, rng);
      log_inplace(t, n, rng);
      for (size_t j = 0; j < n; j++) t[j] *= inv_shape;
      exp_inplace(t, n, rng);
      for (size_t j = 0; j < n; j++) x[i + j] *= t[j];
      i += n;
    }
  }
}

ALWAYS_INLINE void gen_beta(double *x, size_t len, double a, double b, randompack_rng
									 *rng) {
  // X = Ga/(Ga + Gb) ~ Beta(a,b) if Ga ~ Gamma(a,1) and Gb ~ Gamma(b,1)
  gen_gamma(x, len, a, 1.0, rng);
  double gb[128];
  size_t i = 0;
  while (i < len) {
    int n = mini(128, (int)(len - i));
    gen_gamma(gb, n, b, 1.0, rng);
    double hi = 1.0 - DBL_EPSILON;
    for (int j = 0; j < n; j++) { // force the beta-draws to (0, 1)
      double ga = x[i + j];
      double gbj = gb[j];
      double y = ga / (ga + gbj);
      if (y < DBL_MIN) y = DBL_MIN;
      if (y > hi)      y = hi;
      x[i + j] = y;
    }	 
    i += n;
  }
}

ALWAYS_INLINE void gen_t(double *x, size_t len, double nu, randompack_rng *rng) {
  // Z ~ N(0,1), U ~ Gamma(nu/2, 2) independent => X = Z/sqrt(U/nu) ~ t_nu
  double u[128];
  rand_norm(x, len, rng);
  size_t i = 0;
  while (i < len) {
    int n = mini(128, (int)(len - i));
    gen_gamma(u, n, 0.5*nu, 2.0, rng);
    for (int j = 0; j < n; j++)
      x[i + j] /= sqrt(u[j]/nu);
    i += n;
  }
}

ALWAYS_INLINE void gen_f(double *x, size_t len, double nu1, double nu2,
                                 randompack_rng *rng) {
  // F(nu1,nu2) via gamma ratio:
  //   F = (X1*nu2)/(X2*nu1) where X1 ~ Gamma(nu1/2,1), X2 ~ Gamma(nu2/2,1)
  double x2[128];
  gen_gamma(x, len, nu1/2.0, 1.0, rng);
  size_t i = 0;
  while (i < len) {
    int n = mini((int)(len - i), 128);
    gen_gamma(x2, n, nu2/2.0, 1.0, rng);
    for (int j = 0; j < n; j++)
      x[i + j] = (x[i + j]*nu2)/(x2[j]*nu1);
    i += n;
  }
}

ALWAYS_INLINE void gen_chi2(double *x, size_t len, double nu,
  randompack_rng *rng) {
  gen_gamma(x, len, 0.5*nu, 2, rng);
}

ALWAYS_INLINE void gen_gumbel(double *x, size_t len, double mu, double beta,
  randompack_rng *rng) {
  rand_dble(x, len, rng); // x in [0,1)
  log_inplace(x, len, rng);
  for (size_t i = 0; i < len; i++) x[i] = -x[i];
  log_inplace(x, len, rng);
  shift_scale_double_inplace(x, len, mu, -beta, rng);
}

ALWAYS_INLINE void gen_pareto(double *x, size_t len, double xm, double alpha,
                               randompack_rng *rng) {
  rand_exp(x, len, rng);
  for (size_t i = 0; i < len; i++) x[i] = x[i]/alpha;
  exp_inplace(x, len, rng);
  scale_double_inplace(x, len, xm, rng);
}

ALWAYS_INLINE void gen_weibull(double *x, size_t len, double shape,
  double scale, randompack_rng *rng) {
  rand_exp(x, len, rng); // x[i] = E ~ Exp(1) using ziggurat exponential
  if (shape == 1) {
    if (scale != 1) scale_double_inplace(x, len, scale, rng);
    return;
  }
  if (shape == 2) {
    for (size_t i = 0; i < len; i++) x[i] = sqrt(x[i]);
    if (scale != 1) scale_double_inplace(x, len, scale, rng);
    return;
  }
  double inv_shape = 1/shape;
  log_inplace(x, len, rng);
  scale_double_inplace(x, len, inv_shape, rng);
  exp_inplace(x, len, rng);
  if (scale != 1) scale_double_inplace(x, len, scale, rng);
}

ALWAYS_INLINE void gen_skew_normal(double *x, size_t len, double mu,
  double sigma, double alpha, randompack_rng *rng) {
  double denom = sqrt(1 + alpha*alpha);
  double delta = alpha/denom;
  double scale = 1/denom;
  rand_norm(x, len, rng);
  size_t i = 0;
  double u[128];
  while (i < len) {
    int n = mini(128, (int)(len - i));
    rand_norm(u, n, rng);
    for (int j = 0; j < n; j++) {
      double v = delta*fabs(u[j]) + scale*x[i + j];
      x[i + j] = mu + sigma*v;
    }
    i += n;
  }
}

// NOTE 1: There are situations when Sig is indefinite but close to being positive
//         definite, for example due to rounding errors. To remedy this a pivoted
//         Cholesky factorization of Sig is employed, which may return an L matrix
//         which is not lower and has some trailing columns 0, but does satisfy
//         L·L' = Sig, and can thus be used to generate vectors with the right
//         distribution.
// NOTE 2: In case the calling program needs the Cholesky factor of Sig this may
//         be returned in L by letting it be a d × d array instead of 0 in the
//         call.
// NOTE 3: When randompack_mvn is to be called multiple times for the same
//         covariance it is possible to save execution time by reusing the
//         Cholesky factorization of Sig on all calls but the first. Specify Sig
//         and return L on the first call, and let Sig be null and specify L on
//         subsequent calls. If both Sig and L are null, the function exits with
//         ok = 0.
static bool rand_mvn(char transp, double mu[], double Sig[], int d, size_t len,
                     double X[], int ldX, double L[], randompack_rng *rng) {
  int info, rank;
  bool Lalloc = false;
  bool TRAN = (transp == 'T');
  double *work = 0, *A = 0, *Wrk = 0;
  int *piv = 0;
  if (Sig) { // Sig needs to be factorized
    if (!L && X == 0) return true;
    if (!L) {
      if (!ALLOC(L, d*d)) goto fail;
      Lalloc = true;
    }
    else {
      laset("Upper", d, d, 0.0, 0.0, L, d);
    }
    if (!ALLOC(work, 2*d)) goto fail;
    if (!ALLOC(piv, d))    goto fail;
    lacpy("Lower", d, d, Sig, d, L, d);
    potrf("Low", d, L, d, &info);
    if (info > 0) {
      lacpy("Lower", d, d, Sig, d, L, d);
      pstrf("Low", d, L, d, piv, &rank, 1e-14, work, &info); // PSD Chol fact.
      if (!ALLOC(A, rank*d)) goto fail;
      copy(rank*d, L, 1, A, 1);
      laset("All", d, d, 0, 0, L, d);
      // Undo pivoting: original row piv[j] gets row j.
      for (int j = 0; j < d; j++)
        copy(rank, A + j, d, L + piv[j], d);
      laset("All", d, d-rank, 0, 0, L + d*rank, d);
      FREE(A);
    }
    FREE(piv);
    FREE(work);
  }
  rank = LastNonzeroColumn(d, d, L);
  if (len == 0) {
    if (Lalloc) FREE(L);
    return true;
  }
  int nlen = (int)len;
  if (rank == d) {
    if (TRAN) {
      for (size_t i = 0; i < len; i++)
        rand_norm(X + i*ldX, d, rng);
      trmm("Left", "Lower", "NoT", "NotUdia", d, nlen, 1.0, L, d, X, ldX);
    }
    else {
      for (int i = 0; i < d; i++)
        rand_norm(X + i*ldX, len, rng);
      trmm("Right", "Lower", "T", "NotUdia", nlen, d, 1.0, L, d, X, ldX);
    }
  }
  else { // Singular Sig, use workspace for the N(0,1) randoms
    if (!ALLOC(Wrk, len*rank)) goto fail;
    rand_norm(Wrk, len*rank, rng);
    if (TRAN)
      gemm("NoT", "NoT", d, nlen, rank, 1.0, L, d, Wrk, rank, 0.0, X, ldX);
    else
      gemm("NoT", "T", nlen, d, rank, 1.0, Wrk, nlen, L, d, 0.0, X, ldX);
    FREE(Wrk);
  }
  if (mu) {
    size_t i;
    if (TRAN) {
      for (i=0; i<len; i++) axpy(d, 1.0, mu, 1, X + i*ldX, 1);
    }
    else {
      for (i=0; i<len; i++) axpy(d, 1.0, mu, 1, X + i, ldX);
    }
  }
  if (Lalloc) FREE(L);
  return true;
 fail:
  FREE(Wrk);
  FREE(A);
  FREE(piv);
  FREE(work);
  if (Lalloc && L) FREE(L);
  return false;
}
