// -*- C -*-
#include "chacha-portable.inc"
// Big-endian compatibility helpers for ChaCha20 state and buffer layout.
static inline uint32_t load_u32_le(uint8_t *src) {
  uint32_t b0 = src[0];
  uint32_t b1 = src[1];
  uint32_t b2 = src[2];
  uint32_t b3 = src[3];
  return b0 | (b1 << 8) | (b2 << 16) | (b3 << 24);
}

static inline void store_u32_le(uint8_t *dst, uint32_t v) {
  dst[0] = (uint8_t)v;
  dst[1] = (uint8_t)(v >> 8);
  dst[2] = (uint8_t)(v >> 16);
  dst[3] = (uint8_t)(v >> 24);
}

// Big-endian compatibility: byte-swap 64-bit words for draw_byte order.
static inline uint64_t bswap64(uint64_t x) {
  return ((x & 0x00000000000000FFULL) << 56) |
         ((x & 0x000000000000FF00ULL) << 40) |
         ((x & 0x0000000000FF0000ULL) << 24) |
         ((x & 0x00000000FF000000ULL) << 8) |
         ((x & 0x000000FF00000000ULL) >> 8) |
         ((x & 0x0000FF0000000000ULL) >> 24) |
         ((x & 0x00FF000000000000ULL) >> 40) |
         ((x & 0xFF00000000000000ULL) >> 56);
}

static inline ALIGN64 void fill_chacha(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint8_t *dst = (uint8_t *)(buf);
  uint8_t zeros[BUFSIZE*8] = {0};
  size_t nbytes = len*8;
  // Big-endian compatibility: counter stored in little-endian bytes.
  uint32_t counter = load_u32_le(state->u8 + 44);
  uint8_t *key     = &state->u8[0];
  uint8_t *nonce   = &state->u8[32];
  chacha20_xor_stream(dst, zeros, nbytes, key, nonce, counter);
  uint32_t blocks = (uint32_t)(len/8);
  counter += blocks; // count of 64 byte blocks
  store_u32_le(state->u8 + 44, counter);
  // Big-endian compatibility: swap output words for byte-wise draws.
  if (!is_little_endian()) {
    for (size_t i = 0; i < len; i++)
      buf[i] = bswap64(buf[i]);
  }
}
