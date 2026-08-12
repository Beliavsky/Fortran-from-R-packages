# Changelog

## 0.1.0

Initial modern Fortran/FPM translation of expm 1.0-0 computational code.

- Added real/complex Higham-style matrix exponential.
- Added R-style Pade, legacy Taylor, Sidje/RBS Pade, Ward77, and configurable
  Al-Mohy/Higham exponential kernels.
- Added LAPACK matrix balancing and reversal support.
- Added integer matrix powers.
- Added Sidje/EXPOKIT `exp(t*A)v` Krylov action.
- Added Higham SPS and block-enlargement Frechet derivatives.
- Added exact, 1-norm estimated, and Frobenius estimated condition numbers.
- Added Schur-based matrix square root and inverse-scaling matrix logarithm.
- Added eigen and hybrid eigen/Ward methods.
- Added seven regression tests and two examples.
- Retained original source and licensing/provenance material.
