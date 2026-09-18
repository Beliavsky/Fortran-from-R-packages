// -*- C -*-
#include <assert.h>
// Jump functions for the xor-family. The Functions are from https://prng.di.unimi.it/ and
// slightly edited for style, and the constants are computed with two independent programs
// one at https://github.com/nessan/xoshiro, and a Python program in
// github.com//jonasson2/randompack/misc/jumpgeneral.py. They were also compared with the
// precomputed ones on the Italian website for those given there. Nessan's program does
// not handle xorshift128+ but its 2^64 matched the Italian website and the other two were
// obtained with the same Python program as the rest.
//
// Supported jump exponents 32, 64, 96, 128, 192, 253 (the last for the SIMD streams).

static inline void xoshiro256_step(uint64_t s[4]) {
  uint64_t t = s[1] << 17;
  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];
  s[2] ^= t;
  s[3] = ROTL(s[3], 45);
}

static inline void xoshiro256_jump(uint64_t s[4], int p) {
  static const uint64_t
    jump32[4] =
    { 0x58120d583c112f69ULL, 0x7d8d0632bd08e6acULL,
      0x214fafc0fbdbc208ULL, 0x0e055d3520fdb9d7ULL},
    jump64[4] =
    { 0xb13c16e8096f0754ULL, 0xb60d6c5b8c78f106ULL,
      0x34faff184785c20aULL, 0x12e4a2fbfc19bff9ULL},
    jump96[4] =
    { 0x148c356c3114b7a9ULL, 0xcdb45d7def42c317ULL,
      0xb27c05962ea56a13ULL, 0x31eebb6c82a9615fULL},
    jump128[4] =
    { 0x180ec6d33cfd0abaULL, 0xd5a61266f0c9392cULL,
      0xa9582618e03fc9aaULL, 0x39abdc4529b1661cULL},
    jump192[4] =
    { 0x76e15d3efefdcbbfULL, 0xc5004e441c522fb3ULL,
      0x77710069854ee241ULL, 0x39109bb02acbe635ULL},
    jump253[4] =
    { 0xdfca68648b28c5afULL, 0xb56437fb2b753802ULL,
      0xebb82aacdf6ca80dULL, 0xa170e108788db093ULL};
  const uint64_t *jump = 0;
  switch (p) {
  case 32:  jump = jump32;  break;
  case 64:  jump = jump64;  break;
  case 96:  jump = jump96;  break;
  case 128: jump = jump128; break;
  case 192: jump = jump192; break;
  case 253: jump = jump253; break;
  default:  assert(0); return; }
  uint64_t s0, s1, s2, s3;
  s0 = s1 = s2 = s3 = 0;
  for (int i = 0; i < 4; i++) {
    for (int b = 0; b < 64; b++) {
      if (jump[i] & (1ULL << b)) {
        s0 ^= s[0];
        s1 ^= s[1];
        s2 ^= s[2];
        s3 ^= s[3];
      }
      xoshiro256_step(s);
    }
  }
  s[0] = s0;
  s[1] = s1;
  s[2] = s2;
  s[3] = s3;
}

static inline void xorshift128p_jump(uint64_t s[2], int p) {
  static const uint64_t jump32[2] = {0x55fcb2d81da9d3caULL, 0x43b6c9723b2aa348ULL};
  static const uint64_t jump64[2] = {0x8a5cd789635d2dffULL, 0x121fd2155c472f96ULL};
  static const uint64_t jump96[2] = {0xea61c9f1f13962aeULL, 0xa1fe50ef79cfafb2ULL};
  const uint64_t *jump = 0;
  switch (p) {
  case 32:  jump = jump32; break;
  case 64:  jump = jump64; break;
  case 96:  jump = jump96; break;
  default:  assert(0); return; }
  uint64_t s0 = 0, s1 = 0;
  for (int i = 0; i < 2; i++) {
    for (int b = 0; b < 64; b++) {
      if (jump[i] & (1ULL << b)) {
        s0 ^= s[0];
        s1 ^= s[1];
      }
      uint64_t                                                     // next()
        next_s1 = s[0],                                            //
        next_s0 = s[1];                                            //
      s[0] = next_s0;                                              //
      next_s1 ^= next_s1 << 23;                                    //
      s[1] = next_s1 ^ next_s0 ^ (next_s1 >> 18) ^ (next_s0 >> 5); //
    }
  }
  s[0] = s0;
  s[1] = s1;
}

static inline void xoroshiro128pp_jump(uint64_t s[2], int p) {
  static const uint64_t jump32[2] = {0xfcceec21d5c306d9ULL, 0x2e1bcf52f1051044ULL};
  static const uint64_t jump64[2] = {0x2bd7a6a6e99c2ddcULL, 0x0992ccaf6a6fca05ULL};
  static const uint64_t jump96[2] = {0x360fd5f2cf8d5d99ULL, 0x9c6e6877736c46e3ULL};
  const uint64_t *J = 0;
  switch (p) {
  case 32:  J = jump32; break;
  case 64:  J = jump64; break;
  case 96:  J = jump96; break;
  default:  assert(0); return; }
  uint64_t s0 = 0, s1 = 0;
  for (int i = 0; i < 2; i++) {
    for (int b = 0; b < 64; b++) {
      if (J[i] & (1ULL << b)) {
        s0 ^= s[0];
        s1 ^= s[1];
      }
      uint64_t
        next_s0 = s[0],                                     // next()
        next_s1 = s[1];                                     //
      next_s1 ^= next_s0;                                   //
      s[0] = ROTL(next_s0, 49) ^ next_s1 ^ (next_s1 << 21); //
      s[1] = ROTL(next_s1, 28);                             //
    }
  }
  s[0] = s0;
  s[1] = s1;
}
