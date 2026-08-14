# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational translation of tabuSearch 1.2.0.
- Added native binary tabu search with preliminary, intensification,
  diversification, aspiration, tabu-list, repeated-search, and history logic.
- Added computational summary routines.
- Added explicit portable local RNG state.
- Corrected sampled-neighborhood zero-utility, preallocated-history,
  negative-intensification, and diversification tie-handling edge cases.
- Added strict checked tests and a 100-problem exact linear stress test.
