# Porting notes

The translation follows the upstream numerical decomposition rather than the R object hierarchy.

- `grpreg_preprocess` implements centering/scaling, group normalization, rank handling, and group orthogonalization.
- `grpreg_fit_core` translates group-descent and local-coordinate-descent updates for Gaussian, GLM, and Cox paths.
- `grpreg_api` converts fits back to raw predictor scale and supplies prediction, survival diagnostics, likelihoods, residuals, and IC selection.
- `grpreg_cv`, `grpreg_mfdr`, `grpreg_spline`, and `grpreg_data` translate exported package-level computational helpers.

Screening rules are optimization accelerators and are not used by the Fortran path. All groups remain eligible at every coordinate sweep. This makes the implementation simpler and deterministic at the cost of speed on high-dimensional paths.

The local symmetric eigensolver is used only for small within-group Gram matrices and natural-spline constraints. It avoids a system BLAS/LAPACK requirement and keeps Windows FPM builds self-contained.
