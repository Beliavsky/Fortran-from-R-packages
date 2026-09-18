// -*- C -*-
// ChatGPT's translation of Melissa O'Neill's randutils::seed_seq_fe128 original
// https://www.pcg-random.org/ & https://gist.github.com/imneme/540829265469e673d045
// See also:
//   Melissa E. O’Neill, "PCG: A Family of Better Random Number Generators"
//   Technical Report HMC-CS-2014-0905
//   https://www.pcg-random.org/pdf/hmc-cs-2014-0905.pdf
// MIT-license (as the original).

// In C++:
//   randutils::seed_seq_fe128 ss;
//   ss.seed(begin, end);
//   ss.generate(first, last);
//
// In this C version:
//   seed_seq_fe128 ss;
//   seed_seq_fe128_seed(&ss, input, n_input);
//   seed_seq_fe128_generate(&ss, dest, n_words);

#include <stdint.h>
#include <stddef.h>
#define SEED_SEQ_INIT_A      0x43b0d7e5u
#define SEED_SEQ_MULT_A      0x931e8875u
#define SEED_SEQ_INIT_B      0x8b51f9ddu
#define SEED_SEQ_MULT_B      0x58f38dedu
#define SEED_SEQ_MIX_MULT_L  0xca01f9ddu
#define SEED_SEQ_MIX_MULT_R  0x4973f715u
#define SEED_SEQ_XSHIFT      16  // sizeof(uint32_t)*8/2

typedef struct {
  uint32_t mixer[4];
} seed_seq_fe128;

// hash step used in mix_entropy (updates *hash_const)
static uint32_t seed_seq_hash(uint32_t value, uint32_t *hash_const) {
  value ^= *hash_const;
  *hash_const *= SEED_SEQ_MULT_A;
  value *= *hash_const;
  value ^= value >> SEED_SEQ_XSHIFT;
  return value;
}

// mix step used in mix_entropy
static uint32_t seed_seq_mix(uint32_t x, uint32_t y) {
  uint32_t result = SEED_SEQ_MIX_MULT_L*x - SEED_SEQ_MIX_MULT_R*y;
  result ^= result >> SEED_SEQ_XSHIFT;
  return result;
}

// Core entropy mixer: fills mixer[0..3] from input [begin,end)
static void seed_seq_mix_entropy(uint32_t mixer[4], const uint32_t *begin, const uint32_t
                                 *end) {
  uint32_t hash_const = SEED_SEQ_INIT_A;
  const uint32_t *current = begin;
  // Initialize mixer[i] with hashed inputs (or 0)
  for (int i=0; i<4; i++) {
    uint32_t v = 0u;
    if (current != end) {
      v = *current;
      current++;
    }
    mixer[i] = seed_seq_hash(v, &hash_const);
  }
  // Cross-mix all pairs (src,dest) with src != dest
  for (int si=0; si<4; si++) {
    uint32_t src = mixer[si];
    for (int di=0; di<4; di++) {
      if (di == si) continue;
      uint32_t h = seed_seq_hash(src, &hash_const);
      mixer[di] = seed_seq_mix(mixer[di], h);
    }
  }
  // Fold any remaining inputs across all dests
  while (current != end) {
    uint32_t v = *current;
    current++;
    for (int di=0; di<4; di++) {
      uint32_t h = seed_seq_hash(v, &hash_const);
      mixer[di] = seed_seq_mix(mixer[di], h);
    }
  }
}

// Public seed(): Construct 4-word entropy pool from input
static void seed_seq_fe128_seed(seed_seq_fe128 *ss, uint32_t *input, int n_input) {
  seed_seq_mix_entropy(ss->mixer, input, input + n_input);
}

// Public generate(): turn entropy pool into an output sequence
static void seed_seq_fe128_generate(const seed_seq_fe128 *ss, uint32_t *dest, int n_out) {
  int src_index = 0;
  uint32_t hash_const = SEED_SEQ_INIT_B;
  for (int i=0; i<n_out; i++) {
    uint32_t dataval = ss->mixer[src_index];
    src_index++;
    if (src_index == 4) src_index = 0;
    dataval ^= hash_const;
    hash_const *= SEED_SEQ_MULT_B;
    dataval *= hash_const;
    dataval ^= dataval >> SEED_SEQ_XSHIFT;
    dest[i] = dataval;
  }
}
