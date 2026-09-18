// -*- C -*-

#ifndef RANDOMPACK_BUFFER_DRAW_INC
#define RANDOMPACK_BUFFER_DRAW_INC

#include <string.h>

ALWAYS_INLINE void fill_buf_if_empty(randompack_rng *rng) {
  if (UNLIKELY(rng->buf_word >= BUFSIZE)) {
    rng->buf_word = rng->buf_byte = 0;
    rng->fill(rng->buf.u64, BUFSIZE, &rng->state);
  }
}

ALWAYS_INLINE void fill_if_empty(randompack_rng *rng, int *w) {
  if (UNLIKELY(*w >= BUFSIZE)) {
    *w = 0;
    rng->fill(rng->buf.u64, BUFSIZE, &rng->state);
  }
}

ALWAYS_INLINE void fill_if_empty32(randompack_rng *rng, int *w) {
  if (UNLIKELY(*w >= 2*BUFSIZE)) {
    *w = 0;
    rng->fill(rng->buf.u64, BUFSIZE, &rng->state);
  }
}

ALWAYS_INLINE void align16(randompack_rng *rng) { // align to next 2-byte boundary
  if (rng->buf_byte & 1) rng->buf_byte++;
  if (rng->buf_byte == 8) {
    rng->buf_byte = 0;
    rng->buf_word++;
  }
}

ALWAYS_INLINE void align32(randompack_rng *rng) { // align to next 4-byte boundary
  rng->buf_byte = (rng->buf_byte + 3)/4*4;
  if (rng->buf_byte == 8) {
    rng->buf_byte = 0;
    rng->buf_word++;
  }
}

ALWAYS_INLINE void align64(randompack_rng *rng) { // align to beginning or word
  if (rng->buf_byte > 0) {
    rng->buf_byte = 0;
    rng->buf_word++;
  }
}

ALWAYS_INLINE uint8_t draw_byte(randompack_rng *rng) {
  fill_buf_if_empty(rng);
  uint64_t w = rng->buf.u64[rng->buf_word];
  uint8_t b = (uint8_t)(w >> (8*rng->buf_byte));
  rng->buf_byte++;
  if (rng->buf_byte >= 8) {
    rng->buf_byte = 0;
    rng->buf_word++;
  }
  return b;
}

ALWAYS_INLINE uint16_t draw_u16(randompack_rng *rng, int *w) {
  if (UNLIKELY(*w >= 4*BUFSIZE)) {
    *w = rng->buf_word = 0;
    rng->fill(rng->buf.u64, BUFSIZE, &rng->state);
  }
  return rng->buf.u16[(*w)++];
}

ALWAYS_INLINE uint32_t draw_u32(randompack_rng *rng, int *w) {
  if (UNLIKELY(*w >= 2*BUFSIZE)) {
    *w = 0;
    rng->fill(rng->buf.u64, BUFSIZE, &rng->state);
  }
  return rng->buf.u32[(*w)++];
}

ALWAYS_INLINE uint64_t draw_u64(randompack_rng *rng, int *w) {
  if (UNLIKELY(*w >= BUFSIZE)) {
    *w = 0;
    rng->fill(rng->buf.u64, BUFSIZE, &rng->state);
  }
  return rng->buf.u64[(*w)++];
}

ALWAYS_INLINE double draw_u01(randompack_rng *rng, int *w) {
  return U64_TO_DOUBLE(draw_u64(rng, w));
}

ALWAYS_INLINE float draw_u01f(randompack_rng *rng, int *w) {
  return U32_TO_FLOAT(draw_u32(rng, w));
}

static inline int enter_u16_mode(randompack_rng *rng) {
  align16(rng);
  return 4*rng->buf_word + rng->buf_byte/2;
}

static inline int enter_u32_mode(randompack_rng *rng) {
  align32(rng);
  return 2*rng->buf_word + rng->buf_byte/4;
}

static inline int enter_u64_mode(randompack_rng *rng) {
  align64(rng);
  return rng->buf_word;
}

static inline void exit_u16_mode(randompack_rng *rng, int w) {
  rng->buf_word = w/4;
  rng->buf_byte = (w & 3)*2;
}

static inline void exit_u32_mode(randompack_rng *rng, int w) {
  rng->buf_word = w/2;
  rng->buf_byte = (w & 1)*4;
}

static inline void exit_u64_mode(randompack_rng *rng, int w) {
  rng->buf_word = w;
}

static inline void draw_raw_copy(uint8_t *out, size_t nbytes, randompack_rng *rng) {
  size_t i = 0;
  while (i < nbytes) {
    fill_buf_if_empty(rng);
    int pos = 8*rng->buf_word + rng->buf_byte;
    int avail = 8*BUFSIZE - pos;
    size_t rem = nbytes - i;
    int take = avail;
    if (rem < (size_t)avail) take = (int)rem;
    memcpy(out + i, rng->buf.u8 + pos, (size_t)take);
    pos += take;
    rng->buf_word = pos/8;
    rng->buf_byte = pos & 7;
    i += take;
  }
}

static inline void draw_raw_copy32(uint32_t *out, size_t nwords, randompack_rng *rng,
  int *w) {
  size_t i = 0;
  while (i < nwords) {
    fill_if_empty32(rng, w);
    size_t avail = 2*BUFSIZE - (size_t)(*w);
    size_t rem = nwords - i;
    size_t take = rem < avail ? rem : avail;
    memcpy(out + i, rng->buf.u32 + (size_t)(*w), take*sizeof(uint32_t));
    *w += (int)take;
    i += take;
  }
}

#endif
