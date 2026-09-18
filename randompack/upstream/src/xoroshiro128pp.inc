// -*- C -*-

static inline ALIGN64 void fill_xoro128pp(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  uint64_t s0 = state->u64[0];
  uint64_t s1 = state->u64[1];
  for (size_t i = 0; i < len; i++) {
    out[i] = ROTL(s0 + s1, 17) + s0;
    s1 ^= s0;
    s0 = ROTL(s0, 49) ^ s1 ^ (s1 << 21);
    s1 = ROTL(s1, 28);
  }
  state->u64[0] = s0;
  state->u64[1] = s1;
}
