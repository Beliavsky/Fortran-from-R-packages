// -*- C -*-
// This file is extracted from the Random123 official github repository as follows:
//
// 1. Clone https://github.com/DEShawResearch/Random123.git
// 2. cd to Random123/include/Random123, and run the C preprocessor:
//      clang -E -P -DR123_NO_STDINC philox.h | clang-format > philox.inc
// 3. Edit philox.inc by hand to retain only philox4x64_R and what it depends on.
//    Rename the function to philox4x64 and hardwire R = 10 as in the original.
// 4. Add local struct definitions for the counter and key types.
// 5. Change __inline__ to inline and remove __atribute__((always_inline)
// 6. Add _umul128 to allow compilation on windows
//
// The Random123 commit used is 1a2b3c4d5e (February 2025). BSD license, see
// License-Random123.txt.
//
// Kristján Jónasson, jonasson@hi.is, 2025
#ifdef _MSC_VER
#include <intrin.h>
#endif

typedef struct { uint64_t v[2]; } r123array2x64;
typedef struct { uint64_t v[4]; } r123array4x64;

static inline uint64_t mulhilo64(uint64_t a, uint64_t b, uint64_t *hip) {
  uint64_t hi;
  uint64_t lo;
  MUL64_WIDE(a, b, hi, lo);
  *hip = hi;
  return lo;
}

static inline r123array2x64
_philox4x64bumpkey(r123array2x64 key) {
  key.v[0] += (0x9E3779B97F4A7C15ULL);
  key.v[1] += (0xBB67AE8584CAA73BULL);
  return key;
}

static inline r123array4x64
_philox4x64round(r123array4x64 ctr, r123array2x64 key) {
  uint64_t hi0;
  uint64_t hi1;
  uint64_t lo0 = mulhilo64((0xD2E7470EE14C6C93ULL), ctr.v[0], &hi0);
  uint64_t lo1 = mulhilo64((0xCA5A826395121157ULL), ctr.v[2], &hi1);
  r123array4x64 out = {
      {hi1 ^ ctr.v[1] ^ key.v[0], lo1, hi0 ^ ctr.v[3] ^ key.v[1], lo0}};
  return out;
}

typedef r123array4x64 philox4x64_ctr_t;
typedef r123array2x64 philox4x64_key_t;

static inline philox4x64_ctr_t
philox4x64(philox4x64_ctr_t ctr, philox4x64_key_t key) {
  ctr = _philox4x64round(ctr, key);
  key = _philox4x64bumpkey(key); ctr = _philox4x64round(ctr, key);
  key = _philox4x64bumpkey(key); ctr = _philox4x64round(ctr, key);
  key = _philox4x64bumpkey(key); ctr = _philox4x64round(ctr, key);
  key = _philox4x64bumpkey(key); ctr = _philox4x64round(ctr, key);
  key = _philox4x64bumpkey(key); ctr = _philox4x64round(ctr, key);
  key = _philox4x64bumpkey(key); ctr = _philox4x64round(ctr, key);
  key = _philox4x64bumpkey(key); ctr = _philox4x64round(ctr, key);
  key = _philox4x64bumpkey(key); ctr = _philox4x64round(ctr, key);
  key = _philox4x64bumpkey(key); ctr = _philox4x64round(ctr, key);
  return ctr;
}

static inline ALIGN64 void fill_philox(uint64_t *buf, size_t len,
  randompack_state *state) {
  uint64_t *out = buf;
  philox4x64_ctr_t ctr;
  philox4x64_key_t key;
  copy64(ctr.v, state->u64, 4);
  copy64(key.v, state->u64 + 4, 2);
  for (size_t i = 0; i < len; i += 4) {
    philox4x64_ctr_t block = philox4x64(ctr, key);
    copy64(out + i, block.v, 4);
    ctr.v[0]++;
  }
  copy64(state->u64, ctr.v, 4);
}
