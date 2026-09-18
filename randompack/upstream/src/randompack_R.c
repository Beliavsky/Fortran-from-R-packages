// randompack_R.c

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include "randompack.h"

/* ---------- external pointer finalizer ---------- */
static SEXP rng_tag(void){
  return Rf_install("randompack_rng");
}

randompack_rng *randompack_rng_from_R(SEXP ext) {
  if (TYPEOF(ext) != EXTPTRSXP || R_ExternalPtrTag(ext) != rng_tag()) return 0;
  return R_ExternalPtrAddr(ext);
}

static void rng_finalizer(SEXP ext){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (rng){
    randompack_free(rng);
    R_ClearExternalPtr(ext);
  }
}

static uint32_t read_u32_num(SEXP x, R_xlen_t i, char *name){
  if (TYPEOF(x) == INTSXP){
    int v = INTEGER(x)[i];
    if (v == NA_INTEGER)
      Rf_error("%s must not contain NA", name);
    if (v < 0)
      Rf_error("%s entries must be in [0, 2^32-1]", name);
    return (uint32_t)v;
  }
  if (TYPEOF(x) == REALSXP){
    double d = REAL(x)[i];
    if (!R_FINITE(d))
      Rf_error("%s must contain finite values", name);
    if (d < 0 || d > 4294967295.0)
      Rf_error("%s entries must be in [0, 2^32-1]", name);
    uint32_t v = (uint32_t)d;
    if ((double)v != d)
      Rf_error("%s entries must be whole numbers", name);
    return v;
  }
  Rf_error("%s must be a numeric vector", name);
  return 0;
}

/* ---------- constructors / destructors ---------- */
#include <string.h>

SEXP randompack_create_R(SEXP engine, SEXP bitexact_) {
  if (!Rf_isString(engine) || LENGTH(engine) != 1)
    Rf_error("engine must be a single character string");
  const char *name = CHAR(STRING_ELT(engine, 0));
  int bitexact = Rf_asLogical(bitexact_);
  if (bitexact == NA_LOGICAL)
    Rf_error("bitexact must be TRUE or FALSE");
  randompack_rng *rng = randompack_create(name);
  if (!rng) {
    Rf_error("unknown RNG engine '%s' (check spelling)", name);
  }
  char *err = randompack_last_error(rng);
  if (err && err[0]) {
    Rf_error("%s", err);
  }
  if (bitexact) {
    bool ok = randompack_bitexact(rng, true);
    if (!ok) {
      char *msg = randompack_last_error(rng);
      if (!msg) msg = "randompack_bitexact failed";
      Rf_error("%s", msg);
    }
  }
  SEXP ext = R_MakeExternalPtr(rng, rng_tag(), R_NilValue);
  R_RegisterCFinalizerEx(ext, rng_finalizer, TRUE);
  return ext;
}

SEXP randompack_free_R(SEXP ext){
  rng_finalizer(ext);
  return R_NilValue;
}

SEXP randompack_seed_R(SEXP ext, SEXP seed_, SEXP spawn_key_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int seed = Rf_asInteger(seed_);
  if (seed == NA_INTEGER)
    Rf_error("seed must be a non-NA integer");
  int nkey = 0;
  uint32_t *key = 0;
  if (spawn_key_ != R_NilValue) {
    if (TYPEOF(spawn_key_) != INTSXP)
      Rf_error("spawn_key must be an integer vector");
    nkey = LENGTH(spawn_key_);
    if (nkey > 0) {
      int *src = INTEGER(spawn_key_);
      void *mem = R_alloc(nkey, sizeof(uint32_t));
      key = mem;
      for (int i = 0; i < nkey; i++) {
        if (src[i] == NA_INTEGER)
          Rf_error("spawn_key must not contain NA");
        key[i] = src[i];
      }
    }
  }
  bool ok = randompack_seed(seed, key, nkey, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_seed failed";
    Rf_error("%s", msg);
  }
  return R_NilValue;
}

SEXP randompack_randomize_R(SEXP ext){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  bool ok = randompack_randomize(rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_randomize failed";
    Rf_error("%s", msg);
  }
  return R_NilValue;
}

SEXP randompack_duplicate_R(SEXP ext){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  randompack_rng *out = randompack_duplicate(rng);
  if (!out){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_duplicate failed";
    Rf_error("%s", msg);
  }
  char *err = randompack_last_error(out);
  if (err && err[0])
    Rf_error("%s", err);
  SEXP ext_out = R_MakeExternalPtr(out, rng_tag(), R_NilValue);
  R_RegisterCFinalizerEx(ext_out, rng_finalizer, TRUE);
  return ext_out;
}

SEXP randompack_full_mantissa_R(SEXP ext, SEXP enable_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int v = Rf_asLogical(enable_);
  if (v == NA_LOGICAL)
    Rf_error("enable must be TRUE or FALSE");
  bool ok = randompack_full_mantissa(rng, v != 0);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_full_mantissa failed";
    Rf_error("%s", msg);
  }
  return R_NilValue;
}

SEXP randompack_jump_R(SEXP ext, SEXP p_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int p = Rf_asInteger(p_);
  if (p == NA_INTEGER)
    Rf_error("p must be a single integer");
  bool ok = randompack_jump(p, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_jump failed";
    Rf_error("%s", msg);
  }
  return R_NilValue;
}

SEXP randompack_advance_R(SEXP ext, SEXP delta_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  if (!Rf_isNumeric(delta_))
    Rf_error("delta must be a numeric vector");
  R_xlen_t n = XLENGTH(delta_);
  if (n <= 0 || n > 4)
    Rf_error("delta must have length between 1 and 4");
  uint64_t delta[2] = {0, 0};
  for (R_xlen_t i = 0; i < n; i++) {
    uint32_t w = read_u32_num(delta_, i, "delta");
    if (i < 2) {
      if (i & 1) delta[0] |= ((uint64_t)w << 32);
      else delta[0] |= w;
    }
    else {
      if (i & 1) delta[1] |= ((uint64_t)w << 32);
      else delta[1] |= w;
    }
  }
  bool ok = randompack_advance(delta, rng);
  if (!ok) {
    char *err = randompack_last_error(rng);
    Rf_error("%s", err ? err : "randompack_advance failed");
  }
  return R_NilValue;
}


SEXP randompack_set_state_R(SEXP ext, SEXP state_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  if (!Rf_isNumeric(state_))
    Rf_error("state must be a numeric vector");
  R_xlen_t m = XLENGTH(state_);
  if (m <= 0 || (m & 1))
    Rf_error("state must have even length (32-bit words)");
  int nstate = (int)(m/2);
  if ((R_xlen_t)(2*nstate) != m)
    Rf_error("state too long");
  void *mem = R_alloc((size_t)nstate, sizeof(uint64_t));
  uint64_t *u = mem;
  for (int i = 0; i < nstate; i++) {
    uint32_t lo = read_u32_num(state_, 2*i, "state");
    uint32_t hi = read_u32_num(state_, 2*i + 1, "state");
    u[i] = (uint64_t)lo | ((uint64_t)hi << 32);
  }
  bool ok = randompack_set_state(u, nstate, rng);
  if (!ok) {
    char *err = randompack_last_error(rng);
    Rf_error("%s", err ? err : "randompack_set_state failed");
  }
  return R_NilValue;
}

SEXP randompack_pcg64_set_inc_R(SEXP ext, SEXP inc_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  if (!Rf_isNumeric(inc_))
    Rf_error("inc must be a numeric vector");
  R_xlen_t n = XLENGTH(inc_);
  if (n <= 0 || n > 4)
    Rf_error("inc must have length between 1 and 4");
  uint64_t inc[2] = {0, 0};
  for (R_xlen_t i = 0; i < n; i++) {
    uint32_t w = read_u32_num(inc_, i, "inc");
    if (i < 2) {
      if (i & 1) inc[0] |= ((uint64_t)w << 32);
      else inc[0] |= w;
    }
    else {
      if (i & 1) inc[1] |= ((uint64_t)w << 32);
      else inc[1] |= w;
    }
  }
  bool ok = randompack_pcg64_set_inc(inc, rng);
  if (!ok) {
    char *err = randompack_last_error(rng);
    Rf_error("%s", err ? err : "randompack_pcg64_set_inc failed");
  }
  return R_NilValue;
}

SEXP randompack_cwg128_set_weyl_R(SEXP ext, SEXP weyl_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  if (!Rf_isNumeric(weyl_))
    Rf_error("weyl must be a numeric vector");
  R_xlen_t n = XLENGTH(weyl_);
  if (n <= 0 || n > 4)
    Rf_error("weyl must have length between 1 and 4");
  uint64_t weyl[2] = {0, 0};
  for (R_xlen_t i = 0; i < n; i++) {
    uint32_t w = read_u32_num(weyl_, i, "weyl");
    if (i < 2) {
      if (i & 1) weyl[0] |= ((uint64_t)w << 32);
      else weyl[0] |= w;
    }
    else {
      if (i & 1) weyl[1] |= ((uint64_t)w << 32);
      else weyl[1] |= w;
    }
  }
  bool ok = randompack_cwg128_set_weyl(weyl, rng);
  if (!ok) {
    char *err = randompack_last_error(rng);
    Rf_error("%s", err ? err : "randompack_cwg128_set_weyl failed");
  }
  return R_NilValue;
}

SEXP randompack_sfc64_set_abc_R(SEXP ext, SEXP abc_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  if (!Rf_isNumeric(abc_))
    Rf_error("abc must be a numeric vector");
  R_xlen_t n = XLENGTH(abc_);
  if (n <= 0 || n > 6)
    Rf_error("abc must have length between 1 and 6");
  uint64_t abc[3] = {0, 0, 0};
  for (R_xlen_t i = 0; i < n; i++) {
    uint32_t w = read_u32_num(abc_, i, "abc");
    int j = (int)(i/2);
    if (i & 1) abc[j] |= ((uint64_t)w << 32);
    else abc[j] |= w;
  }
  bool ok = randompack_sfc64_set_abc(abc, rng);
  if (!ok) {
    char *err = randompack_last_error(rng);
    Rf_error("%s", err ? err : "randompack_sfc64_set_abc failed");
  }
  return R_NilValue;
}

SEXP randompack_chacha_set_nonce_R(SEXP ext, SEXP nonce_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  if (!Rf_isNumeric(nonce_))
    Rf_error("nonce must be a numeric vector");
  R_xlen_t n = XLENGTH(nonce_);
  if (n <= 0 || n > 3)
    Rf_error("nonce must have length between 1 and 3");
  uint32_t nonce[3] = {0, 0, 0};
  for (R_xlen_t i = 0; i < n; i++)
    nonce[i] = read_u32_num(nonce_, i, "nonce");
  bool ok = randompack_chacha_set_nonce(nonce, rng);
  if (!ok) {
    char *err = randompack_last_error(rng);
    Rf_error("%s", err ? err : "randompack_chacha_set_nonce failed");
  }
  return R_NilValue;
}

SEXP randompack_philox_set_key_R(SEXP ext, SEXP key_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  if (!Rf_isNumeric(key_))
    Rf_error("key must be a numeric vector");
  R_xlen_t key_n = XLENGTH(key_);
  if (key_n <= 0 || key_n > 4)
    Rf_error("key must have length between 1 and 4");
  uint64_t key[2] = {0, 0};
  for (R_xlen_t i = 0; i < key_n; i++) {
    uint32_t w = read_u32_num(key_, i, "key");
    int j = (int)(i/2);
    if (i & 1) key[j] |= ((uint64_t)w << 32);
    else key[j] |= w;
  }
  bool ok = randompack_philox_set_key(key, rng);
  if (!ok) {
    char *err = randompack_last_error(rng);
    Rf_error("%s", err ? err : "randompack_philox_set_key failed");
  }
  return R_NilValue;
}

SEXP randompack_squares_set_key_R(SEXP ext, SEXP key_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  if (!Rf_isNumeric(key_))
    Rf_error("key must be a numeric vector");
  R_xlen_t key_n = XLENGTH(key_);
  if (key_n <= 0 || key_n > 2)
    Rf_error("key must have length between 1 and 2");
  uint64_t key = 0;
  for (R_xlen_t i = 0; i < key_n; i++) {
    uint32_t w = read_u32_num(key_, i, "key");
    if (i & 1) key |= ((uint64_t)w << 32);
    else key |= w;
  }
  bool ok = randompack_squares_set_key(key, rng);
  if (!ok) {
    char *err = randompack_last_error(rng);
    Rf_error("%s", err ? err : "randompack_squares_set_key failed");
  }
  return R_NilValue;
}

/* ---------- draws ---------- */
SEXP randompack_unif_R(SEXP ext, SEXP n_, SEXP a_, SEXP b_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double a = Rf_asReal(a_);
  double b = Rf_asReal(b_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok;
  if (a == 0 && b == 1) ok = randompack_u01(REAL(out), n, rng);
  else ok = randompack_unif(REAL(out), n, a, b, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_unif failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_normal_R(SEXP ext, SEXP n_, SEXP mu_, SEXP sigma_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double mu = Rf_asReal(mu_);
  double sigma = Rf_asReal(sigma_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok;
  if (mu == 0 && sigma == 1) ok = randompack_norm(REAL(out), n, rng);
  else ok = randompack_normal(REAL(out), n, mu, sigma, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_normal failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_skew_normal_R(SEXP ext, SEXP n_, SEXP mu_, SEXP sigma_, SEXP alpha_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double mu = Rf_asReal(mu_);
  double sigma = Rf_asReal(sigma_);
  double alpha = Rf_asReal(alpha_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_skew_normal(REAL(out), n, mu, sigma, alpha, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_skew_normal failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_lognormal_R(SEXP ext, SEXP n_, SEXP mu_, SEXP sigma_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double mu = Rf_asReal(mu_);
  double sigma = Rf_asReal(sigma_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_lognormal(REAL(out), n, mu, sigma, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_lognormal failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_gumbel_R(SEXP ext, SEXP n_, SEXP mu_, SEXP beta_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double mu = Rf_asReal(mu_);
  double beta = Rf_asReal(beta_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_gumbel(REAL(out), n, mu, beta, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_gumbel failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_pareto_R(SEXP ext, SEXP n_, SEXP xm_, SEXP alpha_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double xm = Rf_asReal(xm_);
  double alpha = Rf_asReal(alpha_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_pareto(REAL(out), n, xm, alpha, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_pareto failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_exp_R(SEXP ext, SEXP n_, SEXP scale_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double scale = Rf_asReal(scale_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_exp(REAL(out), n, scale, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_exp failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_gamma_R(SEXP ext, SEXP n_, SEXP shape_, SEXP scale_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double shape = Rf_asReal(shape_);
  double scale = Rf_asReal(scale_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_gamma(REAL(out), n, shape, scale, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_gamma failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_chi2_R(SEXP ext, SEXP n_, SEXP nu_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double nu = Rf_asReal(nu_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_chi2(REAL(out), n, nu, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_chi2 failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_beta_R(SEXP ext, SEXP n_, SEXP a_, SEXP b_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double a = Rf_asReal(a_);
  double b = Rf_asReal(b_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_beta(REAL(out), n, a, b, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_beta failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_t_R(SEXP ext, SEXP n_, SEXP nu_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double nu = Rf_asReal(nu_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_t(REAL(out), n, nu, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_t failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_f_R(SEXP ext, SEXP n_, SEXP nu1_, SEXP nu2_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double nu1 = Rf_asReal(nu1_);
  double nu2 = Rf_asReal(nu2_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_f(REAL(out), n, nu1, nu2, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_f failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_weibull_R(SEXP ext, SEXP n_, SEXP shape_, SEXP scale_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  double shape = Rf_asReal(shape_);
  double scale = Rf_asReal(scale_);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  bool ok = randompack_weibull(REAL(out), n, shape, scale, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_weibull failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

// ---------- discrete ----------
SEXP randompack_int_R(SEXP ext, SEXP n_, SEXP min_, SEXP max_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  int min = Rf_asInteger(min_);
  int max = Rf_asInteger(max_);
  SEXP out = PROTECT(Rf_allocVector(INTSXP, n));
  bool ok = randompack_int(INTEGER(out), n, min, max, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_int failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

SEXP randompack_perm_R(SEXP ext, SEXP n_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  SEXP out = PROTECT(Rf_allocVector(INTSXP, n));
  int *x = INTEGER(out);
  bool ok = randompack_perm(x, n, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_perm failed";
    Rf_error("%s", msg);
  }
  for (int i = 0; i < n; i++) x[i] += 1;
  UNPROTECT(1);
  return out;
}

SEXP randompack_sample_R(SEXP ext, SEXP n_, SEXP k_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  int k = Rf_asInteger(k_);
  if (k < 0) Rf_error("k must be non-negative");
  SEXP out = PROTECT(Rf_allocVector(INTSXP, k));
  int *x = INTEGER(out);
  bool ok = randompack_sample(x, n, k, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_sample failed";
    Rf_error("%s", msg);
  }
  for (int i = 0; i < k; i++) x[i] += 1;
  UNPROTECT(1);
  return out;
}

SEXP randompack_raw_R(SEXP ext, SEXP n_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  SEXP out = PROTECT(Rf_allocVector(RAWSXP, n));
  bool ok = randompack_uint8(RAW(out), n, 0, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_uint8 failed";
    Rf_error("%s", msg);
  }
  UNPROTECT(1);
  return out;
}

// ---------- multivariate ----------
SEXP randompack_mvn_R(SEXP ext, SEXP n_, SEXP Sig_, SEXP mu_){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int n = Rf_asInteger(n_);
  if (n < 0) Rf_error("n must be non-negative");
  if (!Rf_isMatrix(Sig_))
    Rf_error("Sigma must be a matrix");
  SEXP dim = Rf_getAttrib(Sig_, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || LENGTH(dim) != 2)
    Rf_error("Sigma must be a matrix");
  int d = INTEGER(dim)[0];
  int d2 = INTEGER(dim)[1];
  if (d <= 0 || d2 != d)
    Rf_error("Sigma must be a non-empty square matrix");
  if (TYPEOF(Sig_) != REALSXP && TYPEOF(Sig_) != INTSXP)
    Rf_error("Sigma must be numeric");
  int nprotect = 0;
  SEXP Sig = Sig_;
  if (TYPEOF(Sig_) != REALSXP) {
    Sig = PROTECT(Rf_coerceVector(Sig_, REALSXP));
    nprotect++;
  }
  double *Sigv = REAL(Sig);
  double *mu = 0;
  SEXP muv = mu_;
  if (mu_ != R_NilValue) {
    if (TYPEOF(mu_) != REALSXP && TYPEOF(mu_) != INTSXP)
      Rf_error("mu must be numeric");
    if (TYPEOF(mu_) != REALSXP) {
      muv = PROTECT(Rf_coerceVector(mu_, REALSXP));
      nprotect++;
    }
    if (XLENGTH(muv) != d)
      Rf_error("mu length must match Sigma");
    mu = REAL(muv);
  }

  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, n, d));
  nprotect++;
  char transp = 'N';
  bool ok = randompack_mvn(&transp, mu, Sigv, d, (size_t)n, REAL(out), n, 0,
    rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_mvn failed";
    UNPROTECT(nprotect);
    Rf_error("%s", msg);
  }
  UNPROTECT(nprotect);
  return out;
}

/* ---------- serialization ---------- */
SEXP randompack_serialize_R(SEXP ext){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  int len = 0;
  bool ok = randompack_serialize(0, &len, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_serialize failed";
    Rf_error("%s", msg);
  }
  if (len < 0) Rf_error("serialization length is negative");
  SEXP raw = PROTECT(Rf_allocVector(RAWSXP, len));
  if (len > 0){
    ok = randompack_serialize(RAW(raw), &len, rng);
    if (!ok){
      char *msg = randompack_last_error(rng);
      if (!msg) msg = "randompack_serialize failed";
      Rf_error("%s", msg);
    }
  }
  UNPROTECT(1);
  return raw;
}

SEXP randompack_deserialize_R(SEXP ext, SEXP raw){
  randompack_rng *rng = (randompack_rng *)R_ExternalPtrAddr(ext);
  if (!rng) Rf_error("RNG pointer is NULL");
  if (TYPEOF(raw) != RAWSXP)
    Rf_error("raw_state must be a raw vector");
  int len = (int)XLENGTH(raw);
  bool ok = randompack_deserialize(RAW(raw), len, rng);
  if (!ok){
    char *msg = randompack_last_error(rng);
    if (!msg) msg = "randompack_deserialize failed";
    Rf_error("%s", msg);
  }
  return R_NilValue;
}

SEXP randompack_engines_R(void){
  int n = 0;
  int emax = 0;
  int dmax = 0;
  if (!randompack_engines(0, 0, &n, &emax, &dmax))
    Rf_error("randompack_engines failed");
  if (n < 0 || emax <= 0 || dmax <= 0)
    Rf_error("randompack_engines returned invalid sizes");
  char *eng = R_alloc((size_t)n*(size_t)emax, 1);
  char *desc = R_alloc((size_t)n*(size_t)dmax, 1);
  if (!randompack_engines(eng, desc, &n, &emax, &dmax))
    Rf_error("randompack_engines failed");
  SEXP engines = PROTECT(Rf_allocVector(STRSXP, n));
  SEXP descriptions = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; i++) {
    SET_STRING_ELT(engines, i, Rf_mkChar(eng + i*emax));
    SET_STRING_ELT(descriptions, i, Rf_mkChar(desc + i*dmax));
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, engines);
  SET_VECTOR_ELT(out, 1, descriptions);
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("engine"));
  SET_STRING_ELT(names, 1, Rf_mkChar("description"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  SEXP rownames = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(rownames)[0] = NA_INTEGER;
  INTEGER(rownames)[1] = -n;
  Rf_setAttrib(out, R_RowNamesSymbol, rownames);
  SEXP cls = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
  Rf_setAttrib(out, R_ClassSymbol, cls);
  UNPROTECT(6);
  return out;
}
