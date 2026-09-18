// -*- C -*-
#include <assert.h>

static const uint64_t jump32[9] = {
  0x225d447b6009ca81ULL, 0x203038dd5c7c7af4ULL, 0x24431478ee6b2a34ULL,
  0x1c5ba10cf5de7ff4ULL, 0x3da7537aab4d5c21ULL, 0xb9bf984f07b9d6eaULL,
  0x58e772114b65b524ULL, 0x0f535213baabdafaULL, 0x6f0cf98612e25c09ULL
};

static const uint64_t jump64[9] = {
  0xb51c1dfe8c39d3ceULL, 0xc5adf3554a5fff69ULL, 0x0166e1f247d5eccaULL,
  0x1e24936423c97857ULL, 0x9fecd738b2dd0192ULL, 0x10d8c6d28f09f108ULL,
  0x02de9777a3c7eec7ULL, 0x5d30dd7d398eb561ULL, 0xa9ea143361405140ULL
};

static const uint64_t jump96[9] = {
  0x9f1c67142c84c502ULL, 0x024d94e3c4b490e8ULL, 0xe9d460859f0659b6ULL,
  0xd697d9321e8373b1ULL, 0x1164275f61142884ULL, 0xd644d1bd1837c737ULL,
  0xad4191bcf0926c6bULL, 0x2624a1b9ef2c42c0ULL, 0xf671bbcee85222abULL
};

static const uint64_t jump128[9] = {
  0xe5e5397b5c20cb13ULL, 0x0d2dc1dfb1338f79ULL, 0x9c80c974601ad658ULL,
  0x90fde743f57480faULL, 0x144fea51e6afe769ULL, 0xc6b1c11018312892ULL,
  0xc9e664ccabd5195cULL, 0x7640d9f7d294b533ULL, 0x48dece251de941e4ULL
};

static const uint64_t jump192[9] = {
  0xcdd87e6f8ff89dceULL, 0xcf580fe4528cecf6ULL, 0x64ebb69a5327b45fULL,
  0xf254c1e441d1568aULL, 0x4e448efc9f9f97b1ULL, 0x1443a786daa27cfbULL,
  0xa4f6b958be0e7e54ULL, 0xebc258b2f86e436aULL, 0xe8331fcf219a53f9ULL
};

static inline void ranlux_jump(uint64_t s[9], int p) {
  const uint64_t *jump = 0;
  uint64_t z[18];
  switch (p) {
    case 32:  jump = jump32;  break;
    case 64:  jump = jump64;  break;
    case 96:  jump = jump96;  break;
    case 128: jump = jump128; break;
    case 192: jump = jump192; break;
    default:  assert(0); return; }
  mul9x9(z, jump, s);
  mod9x9(s, z);
}
