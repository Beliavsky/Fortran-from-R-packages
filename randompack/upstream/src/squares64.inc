// -*- C -*-

static inline uint64_t squares64(uint64_t ctr, uint64_t key) {
  uint64_t x, y, z, t;

  y = x = ctr*key;
  z = y + key;

  x = x*x + y;  x = (x >> 32) | (x << 32);
  x = x*x + z;  x = (x >> 32) | (x << 32);
  x = x*x + y;  x = (x >> 32) | (x << 32);
  t = x = x*x + z;  x = (x >> 32) | (x << 32);

  return t ^ ((x*x + y) >> 32);
}

#define SQUARES_ROUND(outptr) do {              \
   y = x = ctr*key; z = y + key;                \
   x = x*x + y; x = (x >> 32) | (x << 32);      \
   x = x*x + z; x = (x >> 32) | (x << 32);      \
   x = x*x + y; x = (x >> 32) | (x << 32);      \
   t = x = x*x + z; x = (x >> 32) | (x << 32);  \
   outptr = t ^ ((x*x + y) >> 32);              \
   ctr++;                                       \
 } while(0)

static inline ALIGN64 void fill_squares(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  uint64_t ctr = state->u64[0];
  uint64_t key = state->u64[1];
  uint64_t x, y, z, t;
  for (size_t i = 0; i < len; i += 4) {
    y = x = ctr*key; z = y + key;
    x = x*x + y;  x = (x >> 32) | (x << 32);
    x = x*x + z;  x = (x >> 32) | (x << 32);
    x = x*x + y;  x = (x >> 32) | (x << 32);
    t = x = x*x + z;  x = (x >> 32) | (x << 32);
    out[i] = t ^ ((x*x + y) >> 32);
    ctr++;
    y = x = ctr*key; z = y + key;
    x = x*x + y;  x = (x >> 32) | (x << 32);
    x = x*x + z;  x = (x >> 32) | (x << 32);
    x = x*x + y;  x = (x >> 32) | (x << 32);
    t = x = x*x + z;  x = (x >> 32) | (x << 32);
    out[i+1] = t ^ ((x*x + y) >> 32);
    ctr++;
    y = x = ctr*key; z = y + key;
    x = x*x + y;  x = (x >> 32) | (x << 32);
    x = x*x + z;  x = (x >> 32) | (x << 32);
    x = x*x + y;  x = (x >> 32) | (x << 32);
    t = x = x*x + z;  x = (x >> 32) | (x << 32);
    out[i+2] = t ^ ((x*x + y) >> 32);
    ctr++;
    y = x = ctr*key; z = y + key;
    x = x*x + y;  x = (x >> 32) | (x << 32);
    x = x*x + z;  x = (x >> 32) | (x << 32);
    x = x*x + y;  x = (x >> 32) | (x << 32);
    t = x = x*x + z;  x = (x >> 32) | (x << 32);
    out[i+3] = t ^ ((x*x + y) >> 32);
    ctr++;
  }
  state->u64[0] = ctr;
}
