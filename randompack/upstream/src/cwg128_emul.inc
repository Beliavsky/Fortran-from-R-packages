// -*- C -*-
typedef struct {
  uint64_t lo;
  uint64_t hi;
} rp_u128;

#define U128_SET(x, lo_, hi_) do { \
  (x).lo = (lo_);                  \
  (x).hi = (hi_);                  \
} while (0)

#define U128_OR1(x) do { \
  (x).lo |= 1;           \
} while (0)

#define U128_XOR(out, x, y) do { \
  (out).lo = (x).lo ^ (y).lo;    \
  (out).hi = (x).hi ^ (y).hi;    \
} while (0)

#define U128_XOR64(x, v) do { \
  (x).lo ^= (v);              \
} while (0)

#define U128_ADD(x, y) do {                \
  uint64_t _lo = (x).lo + (y).lo;         \
  uint64_t _cy = (_lo < (x).lo);          \
  (x).lo = _lo;                           \
  (x).hi = (x).hi + (y).hi + _cy;         \
} while (0)

#define U128_SHR1(out, x) do {                    \
  (out).lo = ((x).lo >> 1) | ((x).hi << 63);     \
  (out).hi = (x).hi >> 1;                        \
} while (0)

#define U128_SHR96_LO(x) ((x).hi >> 32)

/* Low 128 bits of x*b, where x and b are 128-bit. */
#define U128_MUL(out, x, b) do {                 \
  uint64_t _hi00, _lo00, _hi01, _lo01;          \
  uint64_t _hi10, _lo10;                        \
  MUL64_WIDE((x).lo, (b).lo, _hi00, _lo00);     \
  MUL64_WIDE((x).lo, (b).hi, _hi01, _lo01);     \
  MUL64_WIDE((x).hi, (b).lo, _hi10, _lo10);     \
  (out).lo = _lo00;                             \
  (out).hi = _hi00 + _lo01 + _lo10;             \
} while (0)

static inline ALIGN64 void fill_cwg128(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  // State layout assumed:
  //   state->u64[0:1] : c0 low, c0 high  (must be odd)
  //   state->u64[2:3] : c1 low, c1 high
  //   state->u64[4:5] : c2 low, c2 high
  //   state->u64[6:7] : c3 low, c3 high
  rp_u128 c0, c1, c2, c3;
  U128_SET(c0, state->u64[0], state->u64[1]);
  U128_SET(c1, state->u64[2], state->u64[3]);
  U128_SET(c2, state->u64[4], state->u64[5]);
  U128_SET(c3, state->u64[6], state->u64[7]);
  for (size_t i=0; i + 1 < len; i += 2) {
    rp_u128 t1, t2;
    U128_ADD(c2, c1);
    U128_OR1(c2);
    U128_SHR1(t1, c1);
    U128_MUL(t2, t1, c2);
    U128_ADD(c3, c0);
    U128_XOR(c1, t2, c3);
    rp_u128 z = c1;
    U128_XOR64(z, U128_SHR96_LO(c2));
    out[i] = z.lo;
    out[i + 1] = z.hi;
  }
  state->u64[0] = c0.lo;
  state->u64[1] = c0.hi;
  state->u64[2] = c1.lo;
  state->u64[3] = c1.hi;
  state->u64[4] = c2.lo;
  state->u64[5] = c2.hi;
  state->u64[6] = c3.lo;
  state->u64[7] = c3.hi;
}
