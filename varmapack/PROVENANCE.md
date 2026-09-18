# Provenance

## Upstream

- R package: `varmapack`
- Upstream version: 0.1.1
- Upstream author/copyright holder: Kristján Jónasson
- Upstream license: MIT
- Upstream project: `https://github.com/jonasson2/varmapack`
- Translation basis: the source archive supplied by the user.

Selected upstream R, C, test, and metadata files are retained below
`upstream/`. The original `inst/THIRD-PARTY-NOTICES` is preserved verbatim.
The upstream bundled SLICOT and LAPACK-derived implementation files themselves
are intentionally not copied into the maintained translation.

## Translation approach

The Fortran model uses the same VARMA/VARMAX equations and array orientation as
the R package:

`x_t = A_1 x_(t-1) + ... + A_p x_(t-p) + e_t + B_1 e_(t-1) + ... + B_q e_(t-q)`

with optional exogenous terms `C_k z_(t-k+1)` and Gaussian innovations with
covariance `Sigma`.

The main numerical mappings are:

- sample autocovariances and covariance-to-correlation conversion: direct
  translations of the upstream computational definitions;
- PSI matrices: direct recursion `Psi_0=I`,
  `Psi_j=B_j+sum_i A_i Psi_(j-i)`;
- AR/MA spectral radii: companion matrices evaluated by `rfortran-linalg`;
- theoretical autocovariances: an augmented state-space covariance satisfying
  the discrete Lyapunov equation `P=F P F^T+Q`;
- orthogonalized IRFs: PSI matrices times a covariance square root;
- simulation: exact stationary Gaussian startup distributions (no burn-in),
  supplied-start conditional Gaussian shocks, and direct VARMA/VARMAX forward
  recursions;
- random draws: sibling `randompack` translation;
- built-in testcases: all 16 named cases plus `random`, `deterministic`, and
  `rho` constructors.

## Material differences from R/C implementation

The maintained Fortran code is a typed numerical API rather than an R6/C-ABI
clone. In particular:

- `model%sim` always provides both `x` and `e` output arrays; callers may ignore
  `e` rather than selecting R's `return_shocks` result shape.
- R conditions/exceptions are represented by integer `info` status values.
- R matrices/arrays are represented directly by Fortran rank-2/rank-3 arrays.
- The stationary covariance setup uses a converged state-space Lyapunov
  fixed-point iteration rather than the upstream VYW/SLICOT selection logic.
  This targets the same stationary covariance but can differ in rounding and
  performance for difficult models.
- The singular-PSD IRF covariance fallback uses a symmetric eigensquare root
  if ordinary Cholesky fails, so a singular covariance may have a different
  factor orientation even though `L L^T` reconstructs the covariance.
- Random number generation is delegated to the sibling Fortran `randompack`.
  Equivalent simulation paths are deterministic for a fixed Fortran RNG state,
  but startup factorization/order is not claimed to be bit-identical to the R/C
  package for every model.
- R6 object behavior, R data frames, printing, and other R presentation
  behavior are intentionally omitted.
