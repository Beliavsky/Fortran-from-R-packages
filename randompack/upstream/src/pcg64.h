// From github.com/fanf2/pcg-dxsm, 79eea6ae8e83, 2025-10-29. See also pcg64_dxsm.inc
// The #defines are from pcg64.def and the struct from pcg.h
#define pcg_ulong_t      __uint128_t
#define pcg_t            pcg64_t

typedef struct pcg {
	pcg_ulong_t state, inc;
} pcg_t;
