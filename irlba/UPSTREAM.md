# Upstream provenance

Source translated from the supplied `irlba` R package source tree:

- Package: `irlba`
- Version: 2.3.7
- Date: 2026-01-26
- Authors/copyright holders listed by upstream:
  - Jim Baglama
  - Lothar Reichel
  - B. W. Lewis
- Upstream package license: GPL-3
- Native C source `src/irlb.c` and `src/utility.c` state:
  Copyright (c) 2016 Bryan W. Lewis, GPL version 3 or later.

The Fortran IRLB kernel is derived from the algorithm and implementation in
those supplied C sources. Higher-level routines are translated from the
supplied R sources `eigen.R`, `prcomp.R`, `ssvd.R`, `svdr.R`, and utilities.

References retained from upstream include:

- J. Baglama and L. Reichel, "Augmented Implicitly Restarted Lanczos
  Bidiagonalization Methods," SIAM Journal on Scientific Computing 27(1),
  2005, 19-42.
- N. Halko, P.-G. Martinsson, and J. Tropp, "Finding Structure with
  Randomness: Stochastic Algorithms for Constructing Approximate Matrix
  Decompositions," 2009.
- H. Shen and J. Z. Huang, "Sparse Principal Component Analysis via
  Regularized Low Rank Matrix Approximation," Journal of Multivariate
  Analysis 99(6), 2008, 1015-1034.
