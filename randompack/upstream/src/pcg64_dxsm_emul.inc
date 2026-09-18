// -*- mode: c -*-

// Emulated 128-bit implementation of PCG64-DXSM

// 128x64 -> 128 via two 64x64 multiplies
static inline void mul128_64(uint64_t a_lo, uint64_t a_hi, uint64_t b,
  uint64_t *out_lo, uint64_t *out_hi) {
  uint64_t hi0;
  uint64_t lo0;
  uint64_t hi1;
  uint64_t lo1;
  MUL64_WIDE(a_lo, b, hi0, lo0);
  MUL64_WIDE(a_hi, b, hi1, lo1);
  *out_lo = lo0;
  *out_hi = hi0 + lo1;
}

static inline ALIGN64 void
fill_pcg64(uint64_t *buf, size_t len, randompack_state *state) {
  uint64_t *out = buf;
  // State layout: (state_lo, state_hi, inc_lo, inc_hi)
  uint64_t next_lo = state->pcg.state_lo;
  uint64_t next_hi = state->pcg.state_hi;
  uint64_t inc_lo = state->pcg.inc_lo;
  uint64_t inc_hi = state->pcg.inc_hi;
  const uint64_t mul = 15750249268501108917ULL;
  for (size_t i = 0; i < len; i++) {
    // LCG step in emulated 128-bit arithmetic
    uint64_t state_lo = next_lo;
    uint64_t state_hi = next_hi;
    uint64_t prod_lo;
    uint64_t prod_hi;
    mul128_64(state_lo, state_hi, mul, &prod_lo, &prod_hi);
    add128(prod_lo, prod_hi, inc_lo, inc_hi, &next_lo, &next_hi);
    uint64_t hi = state_hi; // high 64 bits of state
    uint64_t lo = state_lo | 1;
    hi ^= hi >> 32;
    hi *= mul;
    hi ^= hi >> 48;
    hi *= lo;
    out[i] = hi;
  }
  state->pcg.state_lo = next_lo;
  state->pcg.state_hi = next_hi;
  state->pcg.inc_lo = inc_lo;
  state->pcg.inc_hi = inc_hi;
}
