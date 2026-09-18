// -*- C -*-
static inline ALIGN64 void fill_cwg128(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  uint128_t c0 = state->u128[0]; // odd increment
  uint128_t c1 = state->u128[1];
  uint128_t c2 = state->u128[2];
  uint128_t c3 = state->u128[3];
  for (size_t i=0; i<len; i += 2) {
    c1 = (c1 >> 1) * ((c2 += c1) | 1) ^ (c3 += c0);
    uint128_t z = (c2 >> 96) ^ c1;
    out[i] = (uint64_t)z;
    out[i + 1] = (uint64_t)(z >> 64);
  }
  state->u128[0] = c0;
  state->u128[1] = c1;
  state->u128[2] = c2;
  state->u128[3] = c3;
}
