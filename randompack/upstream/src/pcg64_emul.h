// Emulated PCG64 state layout (no 128-bit integer type available)
#define pcg_t            pcg64_t

typedef struct pcg {
  uint64_t state_lo;
  uint64_t state_hi;
  uint64_t inc_lo;
  uint64_t inc_hi;
} pcg_t;
