// -*- C -*-

#if HAVE_NEON
#include <arm_neon.h>
#define VEC_T uint64x2_t
#define VEC_LOAD(p) vld1q_u64((const uint64_t *)(p))
#define VEC_STORE(p, x) vst1q_u64((uint64_t *)(p), (x))
#define VEC_SET1(x) vdupq_n_u64((x))
#define VEC_ADD(a, b) vaddq_u64((a), (b))
#define VEC_XOR(a, b) veorq_u64((a), (b))
#define VEC_SHL3(a) vshlq_n_u64((a), 3)
#define VEC_SHL17(a) vshlq_n_u64((a), 17)
#define VEC_SHR11(a) vshrq_n_u64((a), 11)
#if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
// Compute x * 9 lane-wise as (x << 3) + x with explicit AArch64 AdvSIMD instructions.
static inline VEC_T vec_mul9_asm(VEC_T x) {
  VEC_T res;
  __asm__(
    "shl %0.2d, %1.2d, #3\n\t"
    "add %0.2d, %0.2d, %1.2d"
    : "=&w"(res)
    : "w"(x)
  );
  return res;
}
#define VEC_MUL9(x) vec_mul9_asm((x))
#else
#define VEC_MUL9(x) VEC_ADD(vshlq_n_u64((x), 3), (x))
#endif
#define VEC_ROTL23(x) vorrq_u64(vshlq_n_u64((x), 23), vshrq_n_u64((x), 41))
#define VEC_ROTL24(x) vorrq_u64(vshlq_n_u64((x), 24), vshrq_n_u64((x), 40))
#define VEC_ROTL45(x) vorrq_u64(vshlq_n_u64((x), 45), vshrq_n_u64((x), 19))
#define HAVE_VEC 1
#else
#define HAVE_VEC 0
#endif
