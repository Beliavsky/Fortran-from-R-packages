// -*- mode: c -*-
static inline ALIGN64 void fill_x128p(uint64_t *buf, size_t len, randompack_state *state) {
  // With parameters (23,18,5), from Vigna (2014, arXiv:1404.0390)
  uint64_t *out = buf;
  size_t n = len;
  uint64_t s0 = state->u64[0];
  uint64_t s1 = state->u64[1];
  for (size_t i = 0; i < n; i++) {
    uint64_t t = s1 ^ (s1 << 23);
    uint64_t result = s0 + t;
    s1 = s0;
    s0 = t ^ s0 ^ (t >> 18) ^ (s0 >> 5);
    out[i] = result;
  }
  state->u64[0] = s0;
  state->u64[1] = s1;
}
