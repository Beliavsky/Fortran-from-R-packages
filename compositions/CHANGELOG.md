# Changelog

## 0.3.0

Final numerical-gap pass.

- Added source-pattern cache construction for `gsiCImpAcompClrExpectation`,
  including missing-order maps, conditional predictor matrices, residual
  covariance/Cholesky factors, normalized-residual Monte Carlo draws, and BDL
  acceptance.
- Added both the source ALR-coordinate Monte Carlo covariance and a corrected
  CLR-projected covariance.
- Defined safe behavior for the upstream no-accepted-draw branch, whose C source
  contains an out-of-bounds access.
- Added empirical and maximum Mahalanobis simulation/calibration, corrected and
  uncorrected outlier decisions, and single-component outlier explanation.
- Added native energy-distance k-sample and multivariate-normal GOF tests with
  permutation/parametric-bootstrap inference.
- Compiled the complete supplied tensorA translation into the main library and
  added a `compositions_tensor` bridge for named high-rank tensor algebra and
  statistics.
- Added the six supplied tensorA regression programs to the combined package test
  suite.
- Expanded strict validation to 19 test programs plus the demo and moved the
  project validation standard to Fortran 2018.

## 0.2.0

Parity expansion focused on the remaining numerical algorithms from the initial
translation.

- Added generalized source-style compositional kriging with arbitrary trends,
  row-specific ALR references, and partially observed compositions.
- Added corrected CLR kriging covariance plus a source-centering compatibility
  option for the upstream buggy `err=TRUE` path.
- Added missing-pattern classification/indexing and conditional ALR covariance
  kernels from the experimental `gsiCImpAcomp*` subsystem.
- Added MAR/BDL conditional imputation with truncated Monte Carlo support.
- Added projection regression with explicit source-compatible and corrected
  modes, plus a completed iterative EM-style imputation extension.
- Added `PBhclust` and exhaustive `PBangprox`; all three `gsi.PrinBal` methods
  now have Fortran counterparts.
- Added dedicated regression tests for missing-data kriging, imputation, and
  principal-balance alternatives.

## 0.1.0

Initial modern Fortran/FPM translation of the main computational layer of
`compositions` 2.0-9.
