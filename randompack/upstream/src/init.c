// init.c

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "randompack.h"

// .Call entry points
SEXP randompack_create_R(SEXP engine, SEXP bitexact);
SEXP randompack_seed_R(SEXP ptr, SEXP seed, SEXP spawn_key);
SEXP randompack_randomize_R(SEXP ptr);
SEXP randompack_duplicate_R(SEXP ptr);
SEXP randompack_full_mantissa_R(SEXP ptr, SEXP enable);
SEXP randompack_jump_R(SEXP ptr, SEXP p);
SEXP randompack_advance_R(SEXP ptr, SEXP delta);
SEXP randompack_set_state_R(SEXP ptr, SEXP state);
SEXP randompack_pcg64_set_inc_R(SEXP ptr, SEXP inc);
SEXP randompack_cwg128_set_weyl_R(SEXP ptr, SEXP weyl);
SEXP randompack_sfc64_set_abc_R(SEXP ptr, SEXP abc);
SEXP randompack_chacha_set_nonce_R(SEXP ptr, SEXP nonce);
SEXP randompack_philox_set_key_R(SEXP ptr, SEXP key);
SEXP randompack_squares_set_key_R(SEXP ptr, SEXP key);
SEXP randompack_unif_R(SEXP ptr, SEXP n, SEXP a, SEXP b);
SEXP randompack_normal_R(SEXP ptr, SEXP n, SEXP mu, SEXP sigma);
SEXP randompack_skew_normal_R(SEXP ptr, SEXP n, SEXP mu, SEXP sigma, SEXP alpha);
SEXP randompack_lognormal_R(SEXP ptr, SEXP n, SEXP mu, SEXP sigma);
SEXP randompack_gumbel_R(SEXP ptr, SEXP n, SEXP mu, SEXP beta);
SEXP randompack_pareto_R(SEXP ptr, SEXP n, SEXP xm, SEXP alpha);
SEXP randompack_exp_R(SEXP ptr, SEXP n, SEXP scale);
SEXP randompack_gamma_R(SEXP ptr, SEXP n, SEXP shape, SEXP scale);
SEXP randompack_chi2_R(SEXP ptr, SEXP n, SEXP nu);
SEXP randompack_beta_R(SEXP ptr, SEXP n, SEXP a, SEXP b);
SEXP randompack_t_R(SEXP ptr, SEXP n, SEXP nu);
SEXP randompack_f_R(SEXP ptr, SEXP n, SEXP nu1, SEXP nu2);
SEXP randompack_weibull_R(SEXP ptr, SEXP n, SEXP shape, SEXP scale);
SEXP randompack_mvn_R(SEXP ptr, SEXP n, SEXP Sig, SEXP mu);
SEXP randompack_int_R(SEXP ptr, SEXP n, SEXP min, SEXP max);
SEXP randompack_perm_R(SEXP ptr, SEXP n);
SEXP randompack_sample_R(SEXP ptr, SEXP n, SEXP k);
SEXP randompack_raw_R(SEXP ptr, SEXP n);
SEXP randompack_serialize_R(SEXP ptr);
SEXP randompack_deserialize_R(SEXP ptr, SEXP raw);
SEXP randompack_engines_R(void);
SEXP randompack_free_R(SEXP ptr);
randompack_rng *randompack_rng_from_R(SEXP ext);

static const R_CallMethodDef CallEntries[] = {
  {"randompack_create_R",      (DL_FUNC)&randompack_create_R,      2},
  {"randompack_seed_R",        (DL_FUNC)&randompack_seed_R,        3},
  {"randompack_randomize_R",   (DL_FUNC)&randompack_randomize_R,   1},
  {"randompack_duplicate_R",  (DL_FUNC)&randompack_duplicate_R,  1},
  {"randompack_full_mantissa_R",(DL_FUNC)&randompack_full_mantissa_R, 2},
  {"randompack_jump_R",       (DL_FUNC)&randompack_jump_R,       2},
  {"randompack_advance_R",    (DL_FUNC)&randompack_advance_R,    2},
  {"randompack_set_state_R",   (DL_FUNC)&randompack_set_state_R,   2},
  {"randompack_pcg64_set_inc_R", (DL_FUNC)&randompack_pcg64_set_inc_R, 2},
  {"randompack_cwg128_set_weyl_R", (DL_FUNC)&randompack_cwg128_set_weyl_R, 2},
  {"randompack_sfc64_set_abc_R", (DL_FUNC)&randompack_sfc64_set_abc_R, 2},
  {"randompack_chacha_set_nonce_R", (DL_FUNC)&randompack_chacha_set_nonce_R, 2},
  {"randompack_philox_set_key_R", (DL_FUNC)&randompack_philox_set_key_R, 2},
  {"randompack_squares_set_key_R", (DL_FUNC)&randompack_squares_set_key_R, 2},
  {"randompack_unif_R",        (DL_FUNC)&randompack_unif_R,        4},
  {"randompack_normal_R",      (DL_FUNC)&randompack_normal_R,      4},
  {"randompack_skew_normal_R", (DL_FUNC)&randompack_skew_normal_R, 5},
  {"randompack_lognormal_R",   (DL_FUNC)&randompack_lognormal_R,   4},
  {"randompack_gumbel_R",      (DL_FUNC)&randompack_gumbel_R,      4},
  {"randompack_pareto_R",      (DL_FUNC)&randompack_pareto_R,      4},
  {"randompack_exp_R",         (DL_FUNC)&randompack_exp_R,         3},
  {"randompack_gamma_R",       (DL_FUNC)&randompack_gamma_R,       4},
  {"randompack_chi2_R",        (DL_FUNC)&randompack_chi2_R,        3},
  {"randompack_beta_R",        (DL_FUNC)&randompack_beta_R,        4},
  {"randompack_t_R",           (DL_FUNC)&randompack_t_R,           3},
  {"randompack_f_R",           (DL_FUNC)&randompack_f_R,           4},
  {"randompack_weibull_R",     (DL_FUNC)&randompack_weibull_R,     4},
  {"randompack_mvn_R",         (DL_FUNC)&randompack_mvn_R,         4},
  {"randompack_int_R",         (DL_FUNC)&randompack_int_R,         4},
  {"randompack_perm_R",        (DL_FUNC)&randompack_perm_R,        2},
  {"randompack_sample_R",      (DL_FUNC)&randompack_sample_R,      3},
  {"randompack_raw_R",         (DL_FUNC)&randompack_raw_R,         2},
  {"randompack_serialize_R",   (DL_FUNC)&randompack_serialize_R,   1},
  {"randompack_deserialize_R", (DL_FUNC)&randompack_deserialize_R, 2},
  {"randompack_engines_R",     (DL_FUNC)&randompack_engines_R,     0},
  {"randompack_free_R",        (DL_FUNC)&randompack_free_R,        1},
  {0, 0, 0}
};

void R_init_randompack(DllInfo *dll){
  R_registerRoutines(dll, 0, CallEntries, 0, 0);
  R_RegisterCCallable("randompack", "randompack_create", (DL_FUNC)&randompack_create);
  R_RegisterCCallable("randompack", "randompack_free", (DL_FUNC)&randompack_free);
  R_RegisterCCallable("randompack", "randompack_last_error",
                     (DL_FUNC)&randompack_last_error);
  R_RegisterCCallable("randompack", "randompack_mvn", (DL_FUNC)&randompack_mvn);
  R_RegisterCCallable("randompack", "randompack_seed", (DL_FUNC)&randompack_seed);
  R_RegisterCCallable("randompack", "randompack_u01", (DL_FUNC)&randompack_u01);
  R_RegisterCCallable("randompack", "randompack_rng_from_R",
                     (DL_FUNC)&randompack_rng_from_R);
  R_useDynamicSymbols(dll, FALSE);
}
