// -*- mode: c; -*-
#include <time.h>

static int get_state_words(randompack_rng *rng) {
  for (int i = 0; i < LEN(rng_table); i++) {
    if (rng_table[i].engine == rng->engine) {
      int nwords = rng_table[i].state_words;
      return nwords > 0 ? nwords : 1;
    }
  }
  return 1;
}

static bool select_engine(const char *s, randompack_rng *rng) {
  // set rng->{engine,next} according to the engine name s
  if (!rng) return false;
  if (!s) {
    rng->engine = X256PPSIMD; // default engine
    rng->fill   = fill_x256ppsimd;
    return true;
  }
  char t[64];
  STRSET(t, s);
  for (int i = 0; t[i]; i++) {
    t[i] = TOLOWER(t[i]);
    if (t[i] == '-') t[i] = '_';
  }
  for (int i = 0; i < LEN(rng_table); i++) {
    if (!strcmp(t, rng_table[i].identifier)) {
      rng->engine = rng_table[i].engine;
      rng->fill   = rng_table[i].fill;
      return true;
    }
  }
  return false;  // unknown engine
}

static inline uint64_t rand_splitmix64(uint64_t *x) {
  uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

static void fill_splitmix64(uint8_t *a, int len, uint64_t *x) {
  for (int k = 0; k < len; k += 8) {
    uint64_t mixed = rand_splitmix64(x);
    memcpy(a + k, &mixed, mini(8, len - k));
  }
}

static void fill_entropy(void *buffer, int len) {
  bool ok = crypto_random(buffer, (size_t) len);
  if (!ok) {  // fallback on esoteric systems
    uint64_t x = clock_nsec();
    x ^= (uintptr_t)buffer;
    fill_splitmix64(buffer, len, &x);
  }
}

static inline void get_stream(uint64_t stream_state[4], int i, randompack_rng *rng) {
  xo256 *st = &rng->state.xo;
  stream_state[0] = st->s0[i];
  stream_state[1] = st->s1[i];
  stream_state[2] = st->s2[i];
  stream_state[3] = st->s3[i];
}

static inline void set_stream(int i, const uint64_t stream_state[4], randompack_rng *rng) {
  xo256 *st = &rng->state.xo;
  st->s0[i] = stream_state[0];
  st->s1[i] = stream_state[1];
  st->s2[i] = stream_state[2];
  st->s3[i] = stream_state[3];
}

static inline void scatter_to_stream0(randompack_rng *rng) {
  // used for SIMD streams
  xo256 *st = &rng->state.xo;
  uint64_t *u64 = rng->state.u64;
  st->s0[0] = u64[0];
  st->s1[0] = u64[1];
  st->s2[0] = u64[2];
  st->s3[0] = u64[3];
}

static inline bool is_x256simd_engine(randompack_rng *rng) {
  return rng->engine == X256PPSIMD || rng->engine == X256SSSIMD;
}

static inline void x256simd_jump(int p, randompack_rng *rng) {
  for (int i=0; i<8; i++) {
    uint64_t stream_state[4];
    get_stream(stream_state, i, rng);
    xoshiro256_jump(stream_state, p);
    set_stream(i, stream_state, rng);
  }
}

// › Add, as a convenience, a function pcg_jump implementing this:
//   void pcg_jump(uint128_t *state, uint128_t m, uint128_t inc, int p) {
//     // advance by 2^p steps
//     for (int i = 0; i < p; i++) {
//       inc = (m + 1) * inc;
//       m *= m;
//     }
//     *state = m * (*state) + inc;
//   }

static inline void setup_sfc64simd_streams(randompack_rng *rng) {
  xo256 *st = &rng->state.xo;
  uint64_t a = rng->state.u64[0];
  uint64_t b = rng->state.u64[1];
  uint64_t c = rng->state.u64[2];
  uint64_t s = rng->state.u64[3];
  uint64_t delta = 1ULL << 61;
  for (int i = 0; i < 8; i++) {
    st->s0[i] = a;
    st->s1[i] = b;
    st->s2[i] = c;
    st->s3[i] = s + (uint64_t)i*delta;
  }
}

static inline bool all_zero_state(uint64_t *state, int n) {
  for (int i = 0; i < n; i++)
    if (state[i] != 0) return false;
  return true;
}

static inline void normalize_state(randompack_rng *rng) {
  // Make state nonzero for xor-family and ranlux, odd for cwg128 & pcg64.
  int nwords = get_state_words(rng);
  switch (rng->engine) {
  case X256PP:
  case X256SS:
  case X256SSSIMD:
  case X256PPSIMD:
  case XORO:
  case X128P:
  case RANLUXPP:
    if (nwords > 0 && all_zero_state(rng->state.u64, nwords))
      rng->state.u64[0] = 1;
    break;
  case CWG128:
#if HAVE128
      rng->state.u128[0] |= (uint128_t)1;
#else
      rng->state.u64[0] |= 1;
#endif
      break;
    case PCG64:
#if HAVE128
      rng->state.pcg.inc |= (pcg_ulong_t)1;
#else
      rng->state.u64[2] |= 1;
#endif
    break;
  default:
    break;
  }
}

static inline void fill_stream_states_with_jump(randompack_rng *rng) {
  for (int i=1; i<8; i++) {
    uint64_t stream_state[4];
    get_stream(stream_state, i-1, rng);
    xoshiro256_jump(stream_state, 253);
    set_stream(i, stream_state, rng);
  }
}

static void rand_init(randompack_rng *rng) {
  rng->buf_word = BUFSIZE;
  rng->buf_byte = 0;
  rng->last_error = 0;
}

static void rand_randomize(randompack_rng *rng) {
  if (rng->engine == CHACHA20) {
    uint8_t tmp[44];
    fill_entropy(tmp, sizeof(tmp));
    for (int i = 0; i < 11; i++) {
      uint32_t v = (uint32_t)tmp[4*i] |
        ((uint32_t)tmp[4*i + 1] << 8) |
        ((uint32_t)tmp[4*i + 2] << 16) |
        ((uint32_t)tmp[4*i + 3] << 24);
      rng->state.u32[i] = v;
    }
    rng->state.u32[11] = 0;
  }
  else {
    int nwords = get_state_words(rng);
    fill_entropy(rng->state.u64, nwords*sizeof(uint64_t));
    normalize_state(rng);
    if (is_x256simd_engine(rng)) {
      scatter_to_stream0(rng);
      fill_stream_states_with_jump(rng);
    }
    else if (rng->engine == SFCSIMD) {
      setup_sfc64simd_streams(rng);
    }
  }
}

#include "seed_seq_fe128.inc"
// Build an extended seed sequence: [TAG0, seed32, nkey+1, TAG1, spawn_key[0..nkey-1]]
// then seed seed_seq_fe128 from it and generate nw output words into w. Returns true on
// success, false on invalid args or allocation failure.
#define TAG0 314159265u  // π * 1e8
#define TAG1 271828182u  // e * 1e8
static bool seed_seq_seed(uint32_t *w, int nw, uint32_t seed32, uint32_t *spawn_key, int
                          nkey) {
  if (nkey > 0 && !spawn_key) return false;
  uint32_t *seq = 0;
  int nseq = nkey + 4;
  if (!ALLOC(seq, nseq)) return false;
  seq[0] = TAG0;
  seq[1] = seed32;
  seq[2] = nkey + 1;
  seq[3] = TAG1;
  if (nkey > 0) copy32(seq + 4, spawn_key, nkey);
  seed_seq_fe128 ss;
  seed_seq_fe128_seed(&ss, seq, nseq);
  seed_seq_fe128_generate(&ss, w, nw);
  FREE(seq);
  return true;
}

static void set_state(uint64_t *state, int nstate, randompack_rng *rng) {  
  rng->buf_word = BUFSIZE;
  rng->buf_byte = 0;
  assert(nstate == get_state_words(rng));
  if (rng->engine == CHACHA20) {
    int k = 0;
    for (int i = 0; i < nstate; i++) {
      uint64_t v = state[i];
      rng->state.u32[k++] = (uint32_t)v;
      rng->state.u32[k++] = (uint32_t)(v >> 32);
    }
  }
  else {
    copy64(rng->state.u64, state, nstate);
    if (is_x256simd_engine(rng)) {
      scatter_to_stream0(rng);
      fill_stream_states_with_jump(rng);
    }
    else if (rng->engine == SFCSIMD) {
      setup_sfc64simd_streams(rng);
    }
  }
  rand_init(rng);
}

static void philox_set_key(uint64_t key[2], randompack_rng *rng) {
  rng->state.u64[4] = key[0];
  rng->state.u64[5] = key[1];
  rand_init(rng);
}

static void squares_set_key(uint64_t key, randompack_rng *rng) {
  rng->state.u64[1] = key;
  rand_init(rng);
}

static void sfc64_skip_18_states(randompack_rng *rng) {
  uint64_t discard[18];
  fill_sfc64(discard, LEN(discard), &rng->state);
}

static void sfc64_set_abc(uint64_t abc[3], randompack_rng *rng) {
  rng->state.u64[0] = abc[0];
  rng->state.u64[1] = abc[1];
  rng->state.u64[2] = abc[2];
  sfc64_skip_18_states(rng);
  rand_init(rng);
}

static void pcg64_set_inc(uint64_t inc[2], randompack_rng *rng) {
  rng->state.u64[2] = inc[0];
  rng->state.u64[3] = inc[1];
  rand_init(rng);
}

static void cwg_skip_96_states(randompack_rng *rng) {
  uint64_t discard[192];
  fill_cwg128(discard, LEN(discard), &rng->state);
}

static void cwg128_set_weyl(uint64_t weyl[2], randompack_rng *rng) {
  rng->state.u64[0] = weyl[0];
  rng->state.u64[1] = weyl[1];
  cwg_skip_96_states(rng); // As recommended by Dziala's report
  rand_init(rng);
}

static void chacha_set_nonce(uint32_t nonce[3], randompack_rng *rng) {
  store_u32_le(rng->state.u8 + 32, nonce[0]);
  store_u32_le(rng->state.u8 + 36, nonce[1]);
  store_u32_le(rng->state.u8 + 40, nonce[2]);
  rand_init(rng);
}
