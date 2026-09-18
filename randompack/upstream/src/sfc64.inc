// -*- C -*-
static inline ALIGN64 void fill_sfc64(uint64_t *buf, size_t len, randompack_state *state){
  uint64_t *out = buf;
  uint64_t a = state->u64[0];
  uint64_t b = state->u64[1];
  uint64_t c = state->u64[2];
  uint64_t counter = state->u64[3];
  for (size_t i = 0; i < len; i += 2) {
    uint64_t tmp0 = a + b + counter++;
    uint64_t a0 = b ^ (b >> 11);
    uint64_t b0 = c + (c << 3);
    uint64_t c0 = (c << 24) | (c >> 40);
    c0 += tmp0;
    uint64_t tmp1 = a0 + b0 + counter++;
    a = b0 ^ (b0 >> 11);
    b = c0 + (c0 << 3);
    c = (c0 << 24) | (c0 >> 40);
    c += tmp1;
    out[i] = tmp0;
    out[i+1] = tmp1;
  }
  state->u64[0] = a;
  state->u64[1] = b;
  state->u64[2] = c;
  state->u64[3] = counter;
}
