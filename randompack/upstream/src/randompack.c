// RandomNumbers — random number generation utilities for VARMASIM
//
// See RandomNumbers.h for further information, including parameter descriptions

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "randompack.h"
#include "randompack_internal.h"
#include "BlasGateway.h"
#include "openlibm.inc"
#include "log_exp.inc"
#include "crypto_random.inc"
#include "engines.inc"
#include "buffer_draw.inc"
#include "scale_inplace.inc"

//============================== INTERNAL TYPES AND TABLES ===============================

typedef struct { // Used in table below
  char *identifier;
  char *description;
  rng_engine engine;
  int state_words;
  engine_fill fill;
} rng_entry;

static rng_entry rng_table[] = {  // x256++simd is default
  {"x256++simd","xoshiro256++, SIMD accelerated (8x4x64)", X256PPSIMD, 4,fill_x256ppsimd},
  {"x256**simd","xoshiro256**, SIMD accelerated (8x4x64)", X256SSSIMD, 4,fill_x256sssimd},
  {"sfc64simd","sfc64, SIMD accelerated (8x4x64)",         SFCSIMD,    4,fill_sfc64simd },
  {"x256++",   "xoshiro256++, Vigna & Blackman, 2019 (4x64)",  X256PP,  4,fill_x256pp   },
  {"x256**",   "xoshiro256**, Vigna & Blackman, 2019 (4x64)",  X256SS,  4,fill_x256ss   },
  {"x128+",    "xorshift128+, Vigna, 2014 (2x64)",             X128P,   2,fill_x128p    },
  {"xoro++",   "xoroshiro128++, Vigna & Blackman, 2016 (2x64)",XORO,    2,fill_xoro128pp},
  {"pcg64",    "PCG64-DXSM, O'Neill, 2014 (4x64)",             PCG64,   4,fill_pcg64    },
  {"sfc64",    "sfc64, Chris Doty-Humphrey, 2013 (4x64)",      SFC64,   4,fill_sfc64    },
  {"squares",  "squares64, Widynski, 2021 (2x64)",             SQUARES, 2,fill_squares  },
  {"philox",   "Philox-4x64, Salmon & Moraes, 2011 (6x64)",    PHILOX,  6,fill_philox   },
  {"cwg128",   "cwg128, Działa, 2022 (8x64)",                  CWG128,  8,fill_cwg128   },
  {"ranlux++", "ranlux++, Sibidanov, 2017 (9x64)",             RANLUXPP,9,fill_ranluxpp },
  {"chacha20", "ChaCha20, Bernstein, 2008 (6x64)",             CHACHA20,6,fill_chacha   },
};
// For x256++simd and x256**simd, state.xo stream 0 (4 words) is seeded or initialized
// directly and then jumped by steps of 2^253 to streams 1..7. For sfc64simd, the base
// state words are replicated to 8 streams with counters s + k*2^61 for k = 0..7.

#include "randutil.inc"
#include "distributions.inc"

//=============================== CREATION AND SETUP ===============================

randompack_rng *randompack_create(const char *engine) {
  randompack_rng *rng;
  // Create engine
  if (!ALLOC(rng, 1)) return 0;
  rng->engine = INVALID;
  rng->cpu_has_avx2 = false;
  rng->cpu_has_avx512 = false;
  if (!select_engine(engine, rng)) {
    rng->last_error = "unknown engine name (spelling error in requested engine)";
    return rng;
  }
#if defined(BUILD_AVX512)
  rng->cpu_has_avx512 = cpu_has_avx512();
#endif
#if defined(BUILD_AVX2)
  rng->cpu_has_avx2 = cpu_has_avx2();
  if (rng->engine == X256PPSIMD && rng->cpu_has_avx2) {
    rng->fill = fill_x256ppsimd_avx2;
  }
  if (rng->engine == X256SSSIMD && rng->cpu_has_avx2) rng->fill = fill_x256sssimd_avx2;
  if (rng->engine == SFCSIMD && rng->cpu_has_avx2) rng->fill = fill_sfc64simd_avx2;
#endif
#if defined(BUILD_AVX512)
  if (rng->engine == X256PPSIMD && rng->cpu_has_avx512) {
    rng->fill = fill_x256ppsimd_avx512;
  }
  if (rng->engine == X256SSSIMD && rng->cpu_has_avx512) rng->fill = fill_x256sssimd_avx512;
  if (rng->engine == SFCSIMD && rng->cpu_has_avx512) rng->fill = fill_sfc64simd_avx512;
#endif
  rand_randomize(rng);
  rand_init(rng);
  return rng;
}

bool randompack_seed(int seed, uint32_t *spawn_key, int nkey, randompack_rng *rng)  {
  if (!rng) return false;
  if (rng->engine == INVALID) {
    rng->last_error = "randompack seed: invalid rng";
    return false;
  }
  uint32_t seed32 = seed;
  // Use Melissa O'Neill's seed sequence
  if (nkey < 0 || (nkey > 0 && !spawn_key)) {
    rng->last_error = "randompack seed: invalid spawn_key arguments";
    return false;
  }
  int nwords = get_state_words(rng);
  int nwords32 = nwords * 2;
  uint32_t *w = 0;
  if (!ALLOC(w, nwords32)) {
    rng->last_error = "randompack seed: allocation failed";
    return false;
  }
  bool ok = seed_seq_seed(w, nwords32, seed32, spawn_key, nkey);
  if (!ok) {
    rng->last_error = "randompack seed: allocation failed";
    FREE(w);
    return false;
  }
  copy32(rng->state.u32, w, nwords32);
  FREE(w);
  if (is_x256simd_engine(rng)) {
    scatter_to_stream0(rng);
    fill_stream_states_with_jump(rng);
  }
  else if (rng->engine == SFCSIMD) {
    setup_sfc64simd_streams(rng);
  }
  rand_init(rng);
  return true;
}

bool randompack_randomize(randompack_rng *rng) {
  if (!rng) return false;
  if (rng->engine == INVALID) {
    rng->last_error = "randompack randomize: invalid rng";
    return false;
  }
  rand_randomize(rng);
  rand_init(rng);
  return true;
}

void randompack_free(randompack_rng *rng) {
  if (!rng) return;
  FREE(rng);
}

randompack_rng *randompack_duplicate(randompack_rng *src) {
  if (!src) return 0;
  randompack_rng *dst;
  if (!ALLOC(dst, 1)) return 0;
  memcpy(dst, src, sizeof(*dst));
  dst->last_error = 0;
  return dst;
}

//============================== ENGINE INFORMATION ==============================

bool randompack_engines(char *engines, char *descriptions, int *nengines,
  int *eng_maxlen, int *desc_maxlen) {
  if (!nengines || !eng_maxlen || !desc_maxlen) return false;
  int n = LEN(rng_table);
  int emax = 1;
  int dmax = 1;
  for (int i = 0; i < n; i++) {
    int elen = (int)strlen(rng_table[i].identifier) + 1;
    int dlen = (int)strlen(rng_table[i].description) + 1;
    emax = maxi(emax, elen);
    dmax = maxi(dmax, dlen);
  }
  *nengines = n;
  *eng_maxlen = emax;
  *desc_maxlen = dmax;
  if (!engines) return true;
  if (!descriptions) return false;
  for (int i = 0; i < n; i++) {
    STRSETN(engines + i*emax, emax, rng_table[i].identifier);
    STRSETN(descriptions + i*dmax, dmax, rng_table[i].description);
  }
  return true;
}

char *randompack_last_error(randompack_rng *rng) {
  if (!rng) return 0;
  return rng->last_error;
}

//============================== STREAM SELECTION ==============================

bool randompack_jump(int p, randompack_rng *rng) {
  bool short_jumps;

  if (!rng) return false;
  if (rng->engine == INVALID) {
    rng->last_error = "randompack_jump: invalid rng";
    return false;
  }
  rng->last_error = 0;
  if (rng->engine != X256PP && rng->engine != X256SS &&
      rng->engine != X256SSSIMD && rng->engine != X256PPSIMD &&
      rng->engine != XORO && rng->engine != X128P &&
      rng->engine != PCG64 &&
      rng->engine != RANLUXPP) {
    rng->last_error =
      "randompack_jump: only supported for pcg64, xor-family engines and ranlux++";
    return false;
  }
  short_jumps = (rng->engine == XORO || rng->engine == X128P);
  if (rng->engine == PCG64) {
    if (p < 0 || p > 127) {
      rng->last_error = "unsupported jump exponent (must be in [0,127])";
      return false;
    }
  }
  else if (short_jumps) {
    if (p != 32 && p != 64 && p != 96) {
      rng->last_error = "unsupported jump exponent (must be 32/64/96)";
      return false;
    }
  }
  else if (p != 32 && p != 64 && p != 96 && p != 128 && p != 192) {
    rng->last_error = "unsupported jump exponent (must be 32/64/96/128/192)";
    return false;
  }
  if (rng->engine == X256PP ||
      rng->engine == X256SS)        xoshiro256_jump     (rng->state.u64, p);
  else if (rng->engine == X256SSSIMD) x256simd_jump     (p, rng);
  else if (rng->engine == XORO)     xoroshiro128pp_jump (rng->state.u64, p);
  else if (rng->engine == X128P)    xorshift128p_jump   (rng->state.u64, p);
  else if (rng->engine == RANLUXPP) ranlux_jump         (rng->state.u64, p);
  else if (rng->engine == X256PPSIMD) x256simd_jump     (p, rng);
  else if (rng->engine == PCG64)    pcg_jump            (p, rng);
  rng->buf_word = BUFSIZE;
  rng->buf_byte = 0;
  return true;
}

bool randompack_advance(uint64_t delta[2], randompack_rng *rng) {
  if (!rng) return false;
  if (rng->engine == INVALID) {
    rng->last_error = "randompack advance: invalid rng";
    return false;
  }
  rng->last_error = 0;
  if (!delta) {
    rng->last_error = "randompack advance: delta is null";
    return false;
  }
  if (rng->engine != PCG64) {
    rng->last_error = "randompack advance: only supported for pcg64";
    return false;
  }
  pcg_advance(delta, rng);
  rng->buf_word = BUFSIZE;
  rng->buf_byte = 0;
  return true;
}

bool randompack_pcg64_set_inc(uint64_t inc[2], randompack_rng *rng) {
  if (!rng) return false;
  rng->last_error = 0;
  if (rng->engine != PCG64) {
    rng->last_error = "randompack pcg64_set_inc: engine is not pcg64";
    return false;
  }
  if ((inc[0] & 1) == 0) {
    rng->last_error = "randompack pcg64_set_inc: increment must be odd";
    return false;
  }
  pcg64_set_inc(inc, rng);
  return true;
}

bool randompack_cwg128_set_weyl(uint64_t weyl[2], randompack_rng *rng) {
  if (!rng) return false;
  rng->last_error = 0;
  if (rng->engine != CWG128) {
    rng->last_error = "randompack cwg128_set_weyl: engine is not cwg128";
    return false;
  }
  if ((weyl[0] & 1) == 0) {
    rng->last_error = "randompack cwg128_set_weyl: Weyl increment must be odd";
    return false;
  }
  cwg128_set_weyl(weyl, rng);
  return true;
}

bool randompack_chacha_set_nonce(uint32_t nonce[3], randompack_rng *rng) {
  if (!rng) return false;
  rng->last_error = 0;
  if (rng->engine != CHACHA20) {
    rng->last_error = "randompack chacha_set_nonce: engine is not chacha20";
    return false;
  }
  chacha_set_nonce(nonce, rng);
  return true;
}

bool randompack_philox_set_key(uint64_t key[2], randompack_rng *rng) {
  if (!rng) return false;
  rng->last_error = 0;
  if (rng->engine != PHILOX) {
    rng->last_error = "randompack philox_set_key: engine is not philox";
    return false;
  }
  philox_set_key(key, rng);
  return true;
}

bool randompack_squares_set_key(uint64_t key, randompack_rng *rng) {
  if (!rng) return false;
  rng->last_error = 0;
  if (rng->engine != SQUARES) {
    rng->last_error = "randompack squares_set_key: engine is not squares64";
    return false;
  }
  squares_set_key(key, rng);
  return true;
}

bool randompack_sfc64_set_abc(uint64_t abc[3], randompack_rng *rng) {
  if (!rng) return false;
  rng->last_error = 0;
  if (rng->engine != SFC64) {
    rng->last_error = "randompack sfc64_set_abc: engine is not sfc64";
    return false;
  }
  sfc64_set_abc(abc, rng);
  return true;
}

//================================ CONFIGURATION =================================

bool randompack_full_mantissa(randompack_rng *rng, bool enable) {
  if (!rng) return false;
  if (rng->engine == INVALID) {
    rng->last_error = "randompack full_mantissa: invalid rng";
    return false;
  }
  rng->last_error = 0;
  rng->usefullmantissa = enable;
  return true;
}

bool randompack_bitexact(randompack_rng *rng, bool enable) {
  if (!rng) return false;
  if (rng->engine == INVALID) {
    rng->last_error = "randompack bitexact: invalid rng";
    return false;
  }
  rng->last_error = 0;
  rng->bitexact = enable;
  return true;
}

bool randompack_set_state(uint64_t state[], int nstate, randompack_rng *rng) {
  if (!rng) return false;
  rng->last_error = 0;
  if (!state || nstate < 0) {
    rng->last_error = "randompack set_state: invalid arguments";
    return false;
  }
  if (rng->engine == INVALID) {
    rng->last_error = "randompack set_state: invalid rng";
    return false;
  }
  int nwords = get_state_words(rng);
  if (nwords <= 0)
    rng->last_error = "randompack set_state: unknown engine";
  else if (nstate != nwords)
    rng->last_error = "randompack set_state: wrong nstate for this engine";
  else if (rng->engine == X256PP || rng->engine == X256SS ||
           rng->engine == X256SSSIMD || rng->engine == X256PPSIMD ||
           rng->engine == XORO || rng->engine == X128P ||
           rng->engine == RANLUXPP) {
    if (all_zero_state(state, nstate))
      rng->last_error = "randompack set_state: all-zero state is invalid";
  }
  else if (rng->engine == CWG128 && (state[0] & 1) == 0)
    rng->last_error = "randompack set_state: cwg128 increment must be odd";
  else if (rng->engine == PCG64 && (state[2] & 1) == 0)
    rng->last_error = "randompack set_state: pcg64 increment must be odd";
  if (rng->last_error) return false;
  set_state(state, nstate, rng);
  return true;
}

//================================ SERIALIZATION =================================

static rng_entry *find_entry(rng_engine e) {
  for (int i = 0; i < LEN(rng_table); i++)
    if (rng_table[i].engine == e) return &rng_table[i];
  return 0;
}

enum { RNG_STATE_WORDS = sizeof(((randompack_rng *)0)->state)/sizeof(uint64_t) };

typedef struct {
  uint32_t version;
  rng_engine engine;
  uint64_t state_u64[RNG_STATE_WORDS];
  uint32_t buf_word;
  uint32_t buf_byte;
  uint32_t reserved0;
  uint32_t reserved1;
  uint64_t buf[BUFSIZE];
} rng_blob;

enum {
  STATE_NEED =
  sizeof(uint32_t)*2
  + sizeof(((rng_blob *)0)->state_u64)
  + sizeof(uint32_t)*4
  + sizeof(((rng_blob *)0)->buf)
};

#include "serializations.inc"

bool randompack_serialize(uint8_t *buf, int *len, randompack_rng *rng) {
  // Returns the complete internal state of rng as an opaque byte buffer
  if (!rng) return false;
  rng->last_error = 0;
  if (!len) {
    rng->last_error = "randompack serialize: len is null";
    return false;
  }
  int need = STATE_NEED;
  if (!buf) { // Report needed buffer size
    *len = need;
    return true;
  }
  if (*len < need) {
    rng->last_error = "randompack serialize: buffer too small";
    return false;
  }
  serialize(buf, *len, rng);
  return true;
}

bool randompack_deserialize(const uint8_t *buf, int len, randompack_rng *rng) {
  // Restores the rng state using a buffer obtained with randompack_serialize
  if (!rng) return false;
  rng->last_error = 0;
  if (!buf || len <= 0) {
    rng->last_error = "randompack deserialize: invalid arguments";
    return false;
  }
  rng_blob blob = {0};
  memcpy(&blob, buf, mini(len, STATE_NEED));
  rng_entry *ent = find_entry(blob.engine);
  int need = STATE_NEED;
  if (blob.version != 1 || !ent || len < need) {
    rng->last_error = "randompack deserialize: corrupt state buffer";
    return false;
  }
  if (rng->engine != INVALID && rng->engine != blob.engine) {
    rng->last_error = "randompack deserialize: engine mismatch";
    return false;
  }
  bool ok = deserialize(&blob, ent, rng);
  if (!ok) {
    rng->last_error = "randompack deserialize: allocation failed";
    return false;
  }
#if defined(BUILD_AVX2)
  if (rng->engine == X256PPSIMD && rng->cpu_has_avx2) {
    rng->fill = fill_x256ppsimd_avx2;
  }
  if (rng->engine == X256SSSIMD && rng->cpu_has_avx2) rng->fill = fill_x256sssimd_avx2;
  if (rng->engine == SFCSIMD && rng->cpu_has_avx2) rng->fill = fill_sfc64simd_avx2;
#endif
#if defined(BUILD_AVX512)
  if (rng->engine == X256PPSIMD && rng->cpu_has_avx512) {
    rng->fill = fill_x256ppsimd_avx512;
  }
  if (rng->engine == X256SSSIMD && rng->cpu_has_avx512) rng->fill = fill_x256sssimd_avx512;
  if (rng->engine == SFCSIMD && rng->cpu_has_avx512) rng->fill = fill_sfc64simd_avx512;
#endif
  return true;
}

//=============================== RAW BITSTREAMS ===============================

bool randompack_raw(void *out, size_t nbytes, randompack_rng *rng) {
  if (!rng) return false;
  if (!out)
    rng->last_error = "invalid arguments to randompack_raw";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;
  draw_raw_copy((uint8_t*)out, nbytes, rng);
  return true;
}

bool randompack_uint8(uint8_t x[], size_t len, uint8_t bound, randompack_rng *rng) {
  if (!rng) return false;
  if (!x)
    rng->last_error = "randompack uint8: invalid arguments";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;  
  if (bound == 0)
    draw_raw_copy((uint8_t*)x, len, rng);
  else
	 rand_uint8(x, len, bound, rng);
  return true;
}

bool randompack_uint16(uint16_t x[], size_t len, uint16_t bound, randompack_rng *rng) {
  if (!rng) return false;
  if (!x)
    rng->last_error = "randompack uint16: invalid arguments";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;  
  if (bound == 0) {
    align16(rng);
    draw_raw_copy((uint8_t*)x, 2*len, rng);
  }
  else
	 rand_uint16(x, len, bound, rng);
  return true;
}

bool randompack_uint32(uint32_t x[], size_t len, uint32_t bound, randompack_rng *rng) {
  if (!rng) return false;
  if (!x)
    rng->last_error = "randompack uint32: invalid arguments";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;  
  if (bound == 0) {
    align32(rng);
    draw_raw_copy((uint8_t*)x, 4*len, rng);
  }
  else
	 rand_uint32(x, len, bound, rng);
  return true;
}

bool randompack_uint64(uint64_t x[], size_t len, uint64_t bound, randompack_rng *rng) {
  if (!rng) return false;
  if (!x)
    rng->last_error = "invalid arguments to randompack_uint64";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;
  if (bound == 0) {
    align64(rng);
    draw_raw_copy((uint8_t*)x, 8*len, rng);
  }
  else
	 rand_uint64(x, len, bound, rng);
  return true;
}

//=========================== DISCRETE DISTRIBUTIONS ============================

bool randompack_int(int x[], size_t len, int m, int n, randompack_rng *rng) {
  if (!rng) return false;
  if (!x)
    rng->last_error = "invalid arguments to randompack_int";
  else if (m > n)
    rng->last_error = "randompack int: m must be <= n";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;
  uint32_t bound = (uint32_t)n - (uint32_t)m + 1u;
  if (bound == 0) {
    align32(rng);
    draw_raw_copy((uint8_t*)x, 4*len, rng);
  }
  else
	 rand_int(x, len, m, bound, rng);
  return true;
}

bool randompack_long_long(long long x[], size_t len, long long m, long long n,
                          randompack_rng *rng) {
  if (!rng) return false;
  if (!x)
    rng->last_error = "invalid arguments to randompack_long_long";
  else if (m > n)
    rng->last_error = "randompack long long: m must be <= n";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;
  uint64_t bound = (uint64_t)n - (uint64_t)m + 1u;
  if (bound == 0) {
    align64(rng);
    draw_raw_copy((uint8_t*)x, 8*len, rng);
  }
  else
    rand_long_long(x, len, m, bound, rng);
  return true;
}

bool randompack_perm(int x[], int len, randompack_rng *rng) {
  if (!rng) return false;
  if (!x)
    rng->last_error = "invalid arguments to randompack_perm";
  else if (len > INT_MAX - 1)
    rng->last_error = "randompack perm: len must be <= 2^31 - 2";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;
  rand_perm(x, len, rng);
  return true;
}

bool randompack_sample(int x[], int len, int k, randompack_rng *rng) {
  if (!rng) return false;
  if (!x || k < 0 || k > len)
    rng->last_error = "invalid arguments to randompack_sample";
  else if (len > INT_MAX)
    rng->last_error = "randompack sample: len must be < 2^31";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;
  rand_sample(x, len, k, rng);
  return true;
}

//========================== CONTINUOUS DISTRIBUTIONS ===========================

bool randompack_u01(double x[], size_t len, randompack_rng *rng) {
  if (!rng) return false;
  if (!x)
    rng->last_error = "invalid arguments to randompack_u01";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;
  rand_dble(x, len, rng);
  return true;
}

bool randompack_unif(double x[], size_t len, double a, double b,
  randompack_rng *rng) {
  double w;
  if (!rng) return false;
  if (!x || a >= b) {
    rng->last_error = "invalid arguments to randompack_unif";
    return false;
  }
  rng->last_error = 0;
  if (a==0 && b==1) {
    rand_dble(x, len, rng); // x in [0,1)
    return true;
  }
#if defined(FP_FAST_FMA) // guarantee output <= b when input < 1
  w = nextafter(b - a, 0.0);
#else 
  w = b - a;
#endif
  if (!rng->bitexact && !rng->usefullmantissa) {
#if defined(BUILD_AVX512)
    if (rng->cpu_has_avx512) {
      rand_unif_avx512(x, len, a, w, rng);
      return true;
    }
#endif
#if defined(BUILD_AVX2)
    if (rng->cpu_has_avx2) {
      rand_unif_avx2(x, len, a, w, rng);
      return true;
    }
#endif
#if HAVE_NEON && defined(__aarch64__)
    rand_unif_neon(x, len, a, w, rng);
    return true;
#endif
  }
  rand_dble(x, len, rng); // x in [0,1)
  shift_scale_double_inplace(x, len, a, w, rng);
  return true;
}

bool randompack_norm(double x[], size_t len, randompack_rng *rng) { // standard normal
  if (!rng) return false;
  if (!x)
    rng->last_error = "invalid arguments to randompack_norm";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;
  rand_norm(x, len, rng);
  return true;
}

bool randompack_normal(double x[], size_t len, double mu, double sigma,
  randompack_rng *rng) { // general normal
  if (!rng) return false;
  if (!x || sigma <= 0) {
    rng->last_error = "invalid arguments to randompack_normal";
    return false;
  }
  rng->last_error = 0;
  rand_norm(x, len, rng);
  if (mu != 0 || sigma != 1) shift_scale_double_inplace(x, len, mu, sigma, rng);
  return true;
}

bool randompack_exp(double x[], size_t len, double scale, randompack_rng *rng) {
  if (!rng) return false;
  if (!x || scale <= 0) {
    rng->last_error = "invalid arguments to randompack_exp";
	 return false;
  }
  rng->last_error = 0;
  rand_exp(x, len, rng);
  if (scale != 1) scale_double_inplace(x, len, scale, rng);
  return true;
}

bool randompack_lognormal(double x[], size_t len, double mu, double sigma,
  randompack_rng *rng) {
  if (!rng) return false;
  if (!x || sigma <= 0) {
    rng->last_error = "invalid arguments to randompack_lognormal";
    return false;
  }
  rng->last_error = 0;
  rand_lognormal(x, len, mu, sigma, rng);
  return true;
}

bool randompack_gamma(double x[], size_t len, double shape, double scale,
  randompack_rng *rng) {
  if (!rng) return false;
  if (!x || shape <= 0 || scale <= 0) {
    rng->last_error = "invalid arguments to randompack_gamma";
    return false;
  }
  rng->last_error = 0;
  gen_gamma(x, len, shape, scale, rng);
  return true;
}

bool randompack_beta(double x[], size_t len, double a, double b,
  randompack_rng *rng) {
  if (!rng) return false;
  if (!x || a <= 0 || b <= 0) {
    rng->last_error = "invalid arguments to randompack_beta";
    return false;
  }
  rng->last_error = 0;
  gen_beta(x, len, a, b, rng);
  return true;
}

bool randompack_chi2(double x[], size_t len, double nu, randompack_rng *rng) {
  if (!rng) return false;
  if (!x || nu <= 0) {
    rng->last_error = "invalid arguments to randompack_chi2";
    return false;
  }
  rng->last_error = 0;
  gen_chi2(x, len, nu, rng);
  return true;
}

bool randompack_t(double x[], size_t len, double nu, randompack_rng *rng) {
  if (!rng) return false;
  if (!x || nu <= 0) {
    rng->last_error = "invalid arguments to randompack_t";
    return false;
  }
  rng->last_error = 0;
  gen_t(x, len, nu, rng);
  return true;
}

bool randompack_f(double x[], size_t len, double nu1, double nu2,
  randompack_rng *rng) {
  if (!rng) return false;
  if (!x || nu1 <= 0 || nu2 <= 0) {
    rng->last_error = "invalid arguments to randompack_f";
    return false;
  }
  rng->last_error = 0;
  gen_f(x, len, nu1, nu2, rng);
  return true;
}

bool randompack_gumbel(double x[], size_t len, double mu, double beta,
  randompack_rng *rng) {
  if (!rng) return false;
  if (!x || beta <= 0) {
    rng->last_error = "invalid arguments to randompack_gumbel";
    return false;
  }
  rng->last_error = 0;
  gen_gumbel(x, len, mu, beta, rng);
  return true;
}

bool randompack_pareto(double x[], size_t len, double xm, double alpha,
  randompack_rng *rng) {
  if (!rng) return false;
  if (!x || xm <= 0 || alpha <= 0) {
    rng->last_error = "invalid arguments to randompack_pareto";
    return false;
  }
  rng->last_error = 0;
  gen_pareto(x, len, xm, alpha, rng);
  return true;
}

bool randompack_weibull(double x[], size_t len, double shape, double scale,
  randompack_rng *rng) {
  if (!rng) return false;
  if (!x || shape <= 0 || scale <= 0) {
    rng->last_error = "invalid arguments to randompack_weibull";
    return false;
  }
  rng->last_error = 0;
  gen_weibull(x, len, shape, scale, rng);
  return true;
}

bool randompack_skew_normal(double x[], size_t len, double mu, double sigma,
  double alpha, randompack_rng *rng) {
  if (!rng) return false;
  if (!x || sigma <= 0) {
    rng->last_error = "invalid arguments to randompack_skew_normal";
    return false;
  }
  rng->last_error = 0;
  gen_skew_normal(x, len, mu, sigma, alpha, rng);
  return true;
}

//============================== MULTIVARIATE NORMAL ==============================

bool randompack_mvn(char *transp, double mu[], double Sig[], int d, size_t n, double X[],
                    int ldX, double L[], randompack_rng *rng) {
  if (!rng) return false;
  if (!Sig && !L)
    rng->last_error = "randompack mvn: either Sig or L must be specified";
  else if (!transp || (transp[0] != 'N' && transp[0] != 'T'))
    rng->last_error = "randompack mvn: transp must begin with 'N' or 'T'";
  else if ((!X && n > 0) || d <= 0 || (ldX <= 0 && X))
    rng->last_error = "randompack mvn: invalid arguments";
  else if (X && n > 0 && (size_t)ldX <
           ((transp[0] == 'N') ? n : (size_t)d))
    rng->last_error = "randompack mvn: invalid arguments";
  else
    rng->last_error = 0;
  if (rng->last_error) return false;

  bool ok = rand_mvn(transp[0], mu, Sig, d, n, X, ldX, L, rng);
  if (!ok) {
    rng->last_error = "randompack mvn: memory allocation failed";
    return false;
  }
  else return true;
}
