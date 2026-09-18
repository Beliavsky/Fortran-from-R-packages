# NOTICE and provenance

This directory is a modern free-form Fortran translation of the computational
core of the R package **gamm4**, version 0.3-0.

Upstream package authors and copyright holders include Simon N. Wood, Fabian
Scheipl, and Arno Schneuwly, as identified by the upstream DESCRIPTION file.
The upstream package is licensed under GPL (>= 2). The original package files
used for this translation are preserved under `upstream/`.

Upstream source package:

- R package: `gamm4`
- Version translated: 0.3-0
- Computational source: `upstream/R/gamm4.r`
- CRAN package description date: 2026-08-03

## Reused sibling Fortran packages

The translation intentionally reuses existing top-level packages in
`Fortran-from-R-packages` instead of copying their source:

- `../mgcv` (`mgcv-fortran`) supplies `smooth_spec_t` and the translated smooth
  construction metadata expected by this bridge.
- `../lme4` supplies `random_term_t`, covariance parameterizations, grouped
  random-design construction, and family/covariance constants.
- `../rfortran-linalg` supplies symmetric eigendecomposition, SPD inversion,
  and Cholesky factorization. Its own pinned LAPACK dependency remains external;
  no BLAS or LAPACK source is copied into this package.

No source from these dependencies is vendored in `gamm4`.

## Translation approach

Upstream `gamm4` converts mgcv smooths to random-effect representations and
fits the resulting model with lme4. The Fortran translation follows the same
computational decomposition for explicit numeric designs:

1. diagonalize each single quadratic smooth penalty;
2. retain its null space as fixed effects and scale its penalized space to an
   identity-penalty random-effect block;
3. combine smooth random effects with ordinary lme4-style grouped random
   effects;
4. profile Gaussian ML/REML or use Laplace/PIRLS for supported non-Gaussian
   families; and
5. reconstruct original smooth coefficients, smoothing parameters, covariance,
   and EDF using the `getVb` marginalization formula from upstream.

The pure numerical `getVb` direct path is translated. The optional upstream
reticulate/Python CHOLMOD acceleration is not needed for numerical semantics
and is intentionally omitted.
