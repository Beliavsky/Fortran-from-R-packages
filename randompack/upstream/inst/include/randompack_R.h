#ifndef RANDOMPACK_R_H
#define RANDOMPACK_R_H

#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include "randompack.h"

typedef randompack_rng *(*randompack_rng_from_R_t)(SEXP rng);
typedef randompack_rng *(*randompack_create_t)(const char *engine);
typedef void (*randompack_free_t)(randompack_rng *rng);
typedef char *(*randompack_last_error_t)(randompack_rng *rng);
typedef bool (*randompack_mvn_t)(char *transp, double mu[], double Sig[], int d,
                                 size_t len, double X[], int ldx, double L[],
                                 randompack_rng *rng);
typedef bool (*randompack_seed_t)(int seed, uint32_t *spawn_key, int n_key,
                                  randompack_rng *rng);
typedef bool (*randompack_u01_t)(double x[], size_t len, randompack_rng *rng);

static inline randompack_rng *randompack_rng_from_R(SEXP rng) {
  randompack_rng_from_R_t get_rng =
    (randompack_rng_from_R_t)R_GetCCallable("randompack", "randompack_rng_from_R");
  return get_rng(rng);
}

static inline randompack_rng *randompack_R_create(const char *engine) {
  static randompack_create_t create = 0;
  if (!create) {
    create = (randompack_create_t)R_GetCCallable("randompack", "randompack_create");
  }
  return create(engine);
}

static inline void randompack_R_free(randompack_rng *rng) {
  static randompack_free_t free_rng = 0;
  if (!free_rng) {
    free_rng = (randompack_free_t)R_GetCCallable("randompack", "randompack_free");
  }
  free_rng(rng);
}

static inline char *randompack_R_last_error(randompack_rng *rng) {
  static randompack_last_error_t last_error = 0;
  if (!last_error) {
    last_error = (randompack_last_error_t)
      R_GetCCallable("randompack", "randompack_last_error");
  }
  return last_error(rng);
}

static inline bool randompack_R_mvn(char *transp, double mu[], double Sig[], int d,
                                    size_t len, double X[], int ldx, double L[],
                                    randompack_rng *rng) {
  static randompack_mvn_t mvn = 0;
  if (!mvn) mvn = (randompack_mvn_t)R_GetCCallable("randompack", "randompack_mvn");
  return mvn(transp, mu, Sig, d, len, X, ldx, L, rng);
}

static inline bool randompack_R_seed(int seed, uint32_t *spawn_key, int n_key,
                                     randompack_rng *rng) {
  static randompack_seed_t seed_rng = 0;
  if (!seed_rng) {
    seed_rng = (randompack_seed_t)R_GetCCallable("randompack", "randompack_seed");
  }
  return seed_rng(seed, spawn_key, n_key, rng);
}

static inline bool randompack_R_u01(double x[], size_t len, randompack_rng *rng) {
  static randompack_u01_t u01 = 0;
  if (!u01) u01 = (randompack_u01_t)R_GetCCallable("randompack", "randompack_u01");
  return u01(x, len, rng);
}

#endif // RANDOMPACK_R_H
