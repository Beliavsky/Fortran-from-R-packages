// -*- C -*-

#ifndef RANDOMPACK_INTERNAL_H
#define RANDOMPACK_INTERNAL_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "randompack_config.h"

#if HAVE128
#include "pcg64.h"
#else
#include "pcg64_emul.h"
#endif

typedef struct randompack_rng randompack_rng;

typedef enum {
  INVALID,
  X128P,
  XORO,
  X256SS,
  X256PP,
  X256SSSIMD,
  RANLUXPP,
  SQUARES,
  SFC64,
  PCG64,
  PHILOX,
  CHACHA20,
  CWG128,
  SFCSIMD,
  X256PPSIMD,
} rng_engine;

typedef struct {
  uint64_t s0[8], s1[8], s2[8], s3[8];
} xo256;

typedef union {
  uint8_t u8[48];   // 6 words, used by chacha20
  uint32_t u32[18]; // 9 words, used when seeding and by chacha20
  uint64_t u64[9];  // used by most engines
  xo256 xo;         // 32 words, all used by 8-lane x256++simd/x256**simd/sfc64simd
  pcg64_t pcg;      // 4 words
  #if HAVE128
  uint128_t u128[4];
  #endif
} randompack_state;

typedef void (*engine_fill)(uint64_t *buf, size_t len, randompack_state *state);

void pcg_advance(uint64_t delta[2], randompack_rng *rng);
void pcg_jump(int p, randompack_rng *rng);

struct randompack_rng {
  randompack_state state;
  rng_engine engine;
  int buf_word;
  int buf_byte;
  char *last_error;
  engine_fill fill;
  bool usefullmantissa;
  bool bitexact;
  bool cpu_has_avx2;
  bool cpu_has_avx512;
  union {
    uint8_t u8[8*BUFSIZE];
    uint16_t u16[4*BUFSIZE];
    uint32_t u32[2*BUFSIZE];
    uint64_t u64[BUFSIZE];
  } buf;
};

#if defined(BUILD_AVX512)
bool cpu_has_avx512(void);
void fill_x256ppsimd_avx512(uint64_t *buf, size_t len, randompack_state *state);
void fill_x256sssimd_avx512(uint64_t *buf, size_t len, randompack_state *state);
void fill_sfc64simd_avx512(uint64_t *buf, size_t len, randompack_state *state);
void rand_dble_avx512(double x[], size_t len, randompack_rng *rng);
void rand_unif_avx512(double x[], size_t len, double a, double b,
  randompack_rng *rng);
void rand_float_avx512(float x[], size_t len, randompack_rng *rng);
void scale_double_avx512(double x[], size_t len, double scale);
void scale_float_avx512(float x[], size_t len, float scale);
void shift_scale_double_avx512(double x[], size_t len, double shift, double scale);
void shift_scale_float_avx512(float x[], size_t len, float shift, float scale);
#endif

#if defined(BUILD_AVX2)
bool cpu_has_avx2(void);
void fill_x256ppsimd_avx2(uint64_t *buf, size_t len, randompack_state *state);
void fill_x256sssimd_avx2(uint64_t *buf, size_t len, randompack_state *state);
void fill_sfc64simd_avx2(uint64_t *buf, size_t len, randompack_state *state);
void rand_dble_avx2(double x[], size_t len, randompack_rng *rng);
void rand_unif_avx2(double x[], size_t len, double a, double b,
  randompack_rng *rng);
void rand_float_avx2(float x[], size_t len, randompack_rng *rng);
void scale_double_avx2(double x[], size_t len, double scale);
void scale_float_avx2(float x[], size_t len, float scale);
void shift_scale_double_avx2(double x[], size_t len, double shift, double scale);
void shift_scale_float_avx2(float x[], size_t len, float shift, float scale);
#if defined(RANDOMPACK_TEST_HOOKS)
int randompack_avx2_used(void);
void randompack_avx2_reset(void);
#endif
#endif

#endif
