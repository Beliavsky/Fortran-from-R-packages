// -*- C -*-
#include "randompack_internal.h"

// Advance the PCG64 state by delta. Copied from NumPy
#if HAVE128
void pcg_advance(uint64_t delta[2], randompack_rng *rng) {
  uint128_t state = rng->state.pcg.state;
  uint128_t m = (uint128_t)15750249268501108917ULL;
  uint128_t inc = rng->state.pcg.inc;
  uint128_t d = ((uint128_t)delta[1] << 64) | (uint128_t)delta[0];
  uint128_t m_acc = 1;
  uint128_t inc_acc = 0;
  while (d > 0) {
    if (d & 1) {
      m_acc *= m;
      inc_acc = inc_acc*m + inc;
    }
    inc = (m + 1)*inc;
    m *= m;
    d >>= 1;
  }
  rng->state.pcg.state = m_acc*state + inc_acc;
}

void pcg_jump(int p, randompack_rng *rng) {
  uint128_t state = rng->state.pcg.state;
  uint128_t m = (uint128_t)15750249268501108917ULL;
  uint128_t inc = rng->state.pcg.inc;
  for (int i = 0; i < p; i++) {
    inc = (m + 1)*inc;
    m *= m;
  }
  rng->state.pcg.state = m*state + inc;
}

#else
static inline void mul128_mod(uint64_t a_lo, uint64_t a_hi, uint64_t b_lo,
  uint64_t b_hi, uint64_t *out_lo, uint64_t *out_hi) {
  uint64_t ll_hi, ll_lo, discard_hi, lh_lo, hl_lo;
  MUL64_WIDE(a_lo, b_lo, ll_hi, ll_lo);
  MUL64_WIDE(a_lo, b_hi, discard_hi, lh_lo);
  MUL64_WIDE(a_hi, b_lo, discard_hi, hl_lo);
  *out_lo = ll_lo;
  *out_hi = ll_hi + lh_lo + hl_lo;
}

void pcg_advance(uint64_t delta[2], randompack_rng *rng) {
  uint64_t state_lo = rng->state.pcg.state_lo;
  uint64_t state_hi = rng->state.pcg.state_hi;
  uint64_t m_lo = 15750249268501108917ULL;
  uint64_t m_hi = 0;
  uint64_t inc_lo = rng->state.pcg.inc_lo;
  uint64_t inc_hi = rng->state.pcg.inc_hi;
  uint64_t d_lo = delta[0];
  uint64_t d_hi = delta[1];
  uint64_t m_acc_lo = 1;
  uint64_t m_acc_hi = 0;
  uint64_t inc_acc_lo = 0;
  uint64_t inc_acc_hi = 0;
  while (d_lo != 0 || d_hi != 0) {
    if (d_lo & 1) {
      uint64_t lo, hi;
      mul128_mod(m_acc_lo, m_acc_hi, m_lo, m_hi, &lo, &hi);
      m_acc_lo = lo;
      m_acc_hi = hi;
      mul128_mod(inc_acc_lo, inc_acc_hi, m_lo, m_hi, &lo, &hi);
      add128(lo, hi, inc_lo, inc_hi, &inc_acc_lo, &inc_acc_hi);
    }
    {
      uint64_t mp1_lo, mp1_hi, lo, hi;
      add128(m_lo, m_hi, 1, 0, &mp1_lo, &mp1_hi);
      mul128_mod(mp1_lo, mp1_hi, inc_lo, inc_hi, &lo, &hi);
      inc_lo = lo;
      inc_hi = hi;
      mul128_mod(m_lo, m_hi, m_lo, m_hi, &lo, &hi);
      m_lo = lo;
      m_hi = hi;
    }
    d_lo = (d_lo >> 1) | (d_hi << 63);
    d_hi >>= 1;
  }
  mul128_mod(m_acc_lo, m_acc_hi, state_lo, state_hi, &state_lo, &state_hi);
  add128(state_lo, state_hi, inc_acc_lo, inc_acc_hi, &state_lo, &state_hi);
  rng->state.pcg.state_lo = state_lo;
  rng->state.pcg.state_hi = state_hi;
}

void pcg_jump(int p, randompack_rng *rng) {
  uint64_t state_lo = rng->state.pcg.state_lo;
  uint64_t state_hi = rng->state.pcg.state_hi;
  uint64_t m_lo = 15750249268501108917ULL;
  uint64_t m_hi = 0;
  uint64_t inc_lo = rng->state.pcg.inc_lo;
  uint64_t inc_hi = rng->state.pcg.inc_hi;
  for (int i = 0; i < p; i++) {
    uint64_t mp1_lo, mp1_hi, lo, hi;
    add128(m_lo, m_hi, 1, 0, &mp1_lo, &mp1_hi);
    mul128_mod(mp1_lo, mp1_hi, inc_lo, inc_hi, &lo, &hi);
    inc_lo = lo;
    inc_hi = hi;
    mul128_mod(m_lo, m_hi, m_lo, m_hi, &lo, &hi);
    m_lo = lo;
    m_hi = hi;
  }
  mul128_mod(m_lo, m_hi, state_lo, state_hi, &state_lo, &state_hi);
  add128(state_lo, state_hi, inc_lo, inc_hi, &state_lo, &state_hi);
  rng->state.pcg.state_lo = state_lo;
  rng->state.pcg.state_hi = state_hi;
}
#endif
