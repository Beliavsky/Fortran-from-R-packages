# Porting notes

## Translation policy

The aim is to preserve the computational content of NMOF 2.12-0 while replacing R-specific dynamic objects and external package dispatch with explicit modern Fortran APIs.

The port uses:

- Fortran 2018
- `real(dp)` with `dp = real64`
- `implicit none`
- allocatable arrays and typed results
- abstract interfaces for user callbacks
- explicit status codes rather than R conditions
- a portable explicit RNG state
- BLAS/LAPACK for dense linear algebra

## R-to-Fortran naming

R camel-case names are translated to descriptive snake-case names. Examples:

- `DEopt` -> `de_opt`
- `PSopt` -> `ps_opt`
- `GAopt` -> `ga_opt`
- `vanillaOptionEuropean` -> `vanilla_option_european`
- `mvPortfolio` -> `mean_variance_portfolio`
- `randomReturns` -> `random_returns`
- `repairMatrix` -> `repair_matrix`
- `callHestoncf` -> `call_heston_cf`

The umbrella module `nmof` re-exports the public interfaces of all component modules.

## Numerical substitutions

### Random-number generation

R's global RNG is replaced with a portable explicit generator in `nmof_rng`. Fixed seeds are reproducible within this Fortran implementation, but streams do not match R's `runif` or `rnorm` bit-for-bit.

### Quadratic and linear optimization

The R package can delegate to `quadprog` and `Rglpk`. This port is self-contained apart from BLAS/LAPACK:

- convex quadratic portfolio problems use an active-set QP solver;
- equality, inequality, budget, box, and group constraints use explicit projection and active-set routines;
- CVaR and MAD use deterministic projected subgradient optimization;
- ERC uses deterministic coordinate descent rather than the package's stochastic local-search construction.

The mathematical objectives and constraints are retained, but iteration paths and final rounding may differ from external R solvers.

### Fourier integration

`callCF` and related routines use transformed Gauss-Legendre quadrature rather than R's `integrate`. The characteristic-function formulation is retained.

### Root solving

Yield and implied-volatility solving use bracketed modern Fortran routines where appropriate. The yield-to-maturity Newton convention from NMOF is retained, and `yield_to_maturity_curve` handles vector offsets explicitly.

### Simulation callbacks

The R function `mc(paths, payoff, ...)` only calls the supplied payoff. In Fortran, users call their payoff procedure directly, so no numerical wrapper is necessary.

## Source behavior retained

- Option functions take variance rather than volatility.
- Dividend amounts and continuous dividend yields are mutually exclusive.
- The original Differential Evolution cyclic population shifts are preserved.
- PBO uses average ranks for ties; a numerical tolerance is used for floating-point equality.
- The original moving-average initial-value convention is preserved.
- The original Nelson-Siegel and Nelson-Siegel-Svensson parameter ordering is preserved.

## Corrected implementation defects

The translation corrected implementation defects encountered during validation rather than preserving undefined behavior:

- the bracketed root solver was replaced after the initial implementation failed implied-volatility recovery;
- the PBO combination iterator was corrected for Fortran's non-short-circuit logical evaluation;
- average tie ranks in PBO are calculated explicitly;
- covariance-matrix repair was rewritten to reconstruct the matrix correctly after eigenvalue clipping;
- antithetic and exact-moment random-return allocation paths were corrected;
- projection calls that aliased input and output arrays were separated;
- unavailable American-option Greeks are initialized to IEEE NaN instead of uninitialized storage.

These are translation/integration corrections, not claims that the corresponding current R routines necessarily fail in the same way.

## Features intentionally not translated

The following are not computational algorithms or depend primarily on the R runtime:

- `French`, `Ritter`, and `Shiller`: network downloads and external file parsing;
- `showChapterNames` and `showExample`: package example discovery;
- S3 `print` and `plot` methods;
- `LS.info`, `SA.info`, and `TA.info`: R parent-frame introspection. Fortran callbacks receive iteration counters directly;
- parallel-cluster orchestration, progress bars, and console flushing;
- bundled R data objects and `ts`, list, data-frame, and formula metadata;
- `checkList`, `makeInteger`, `anyNA`, and `mcList`: R argument/list validation infrastructure;
- `mRU` and `mRN`: replaced by the general RNG module;
- `repair1c`: replaced by general bound projection;
- `due`: incorporated into option dividend processing;
- `maxMinR`: the source contains only commented-out unfinished branches and returns an undefined object, so there is no executable algorithm to translate.

## Original source

The unmodified supplied package is retained under `original/`. It is the authoritative provenance reference for names, formulas, defaults, documentation, examples, and licensing.
