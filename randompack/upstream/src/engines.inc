// -*- mode: c; -*-
// RNG engine implementations
#include "xorshift128p.inc"
#include "xoshiro256pp.inc"
#include "xoshiro256ss.inc"
#include "xoroshiro128pp.inc"
#include "simd_macros.inc"
#include "x256pp_neon_scalar.inc"
#include "x256ss_neon_scalar.inc"
#include "sfc64_neon_scalar.inc"
#include "chachacha.inc"
#include "squares64.inc"
#include "sfc64.inc"
#include "jumps.inc"
#include "ranluxpp.inc"
#include "ranlux_jump.inc"

#if HAVE128
#include "pcg64_dxsm.inc"
#include "cwg128.inc"
#else
#include "pcg64_dxsm_emul.inc"
#include "cwg128_emul.inc"
#endif

#include "philox4x64.inc"
