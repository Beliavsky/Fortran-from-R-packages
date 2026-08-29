# Notices and provenance

This directory is a modern Fortran/FPM adaptation of the computational code in
**KFAS 1.6.0**, the R package *Kalman Filter and Smoother for Exponential Family
State Space Models*.

## Upstream

- Upstream author and maintainer: Jouni Helske
- Upstream project: https://github.com/helske/KFAS
- Upstream version translated: 1.6.0
- CRAN publication date in the supplied `DESCRIPTION`: 2025-05-26
- Upstream license declaration: `GPL (>= 2)`
- Upstream ORCID: 0000-0001-7130-793X

The original package metadata, citation record, change log, and SHA-256 hashes of
the supplied computational source are retained under `provenance/`.

## What was retained and adapted

All 30 upstream `src/*.f90` computational source units are retained as the core
of this package. They implement exact diffuse and ordinary Gaussian Kalman
filtering, smoothing, simulation smoothing, non-Gaussian approximation and
importance sampling, weighted simulation moments, state/signal transformations,
and supporting linear-algebra routines.

A modern module, `src/kfas.f90`, provides typed Fortran-facing model, filter,
smoother, and non-Gaussian approximation results. It exposes Gaussian
filtering/smoothing/log-likelihood, automatic correlated-`H` LDL transformation,
deterministic non-Gaussian approximation/log-likelihood, the upstream
`initTheta` calculation, AR transformation, LDL decomposition, and weighted
moments. The specialized simulation, importance-sampling, and additional
non-Gaussian kernels remain part of the library without reimplementing their
algorithms.

## R/Rmath bridge replacement

The upstream `src/cdistwrap.c` file only wrapped Rmath density functions and
`R_FINITE`. It is not copied into this package. `src/rmath_replacements.f90`
implements the same required log-density operations directly in Fortran using
`log`, `log_gamma`, and `ieee_is_finite`. No R headers, R runtime, or Rmath
library is required.

The R registration layer (`init.c`), R object constructors, formula parsing,
S3 methods, plotting, printing, and other R-specific interface code are not
included.

## Shared-dependency review

Before implementation, the top level of
https://github.com/Beliavsky/Fortran-from-R-packages was checked for an existing
KFAS translation and relevant shared numerical modules. No top-level KFAS
translation was present. `rfortran-core` and `rfortran-linalg` were reviewed.
KFAS's retained upstream kernels use the low-level BLAS/LAPACK ABI directly;
the focused shared modules do not provide a drop-in ABI for these calls. The
package therefore declares the MIT-licensed `fortran-lapack` FPM dependency and
uses a small compatibility layer that delegates the retained entry points to
that dependency. It does not copy or vendor BLAS, LAPACK, `r.f90`, `r_mod.f90`,
ARPACK, or any translated R-package dependency.

The ordinary and state-saving Kalman-filter entry points originally contained
parallel implementations. They now share the state-saving implementation, with
compatibility wrappers retaining the original external procedure interfaces.
Maintained numerical declarations and constants use the common `dp = real64`
kind from `kfas_kinds` rather than separate `double precision` declarations and
`d0` literals. This is a source-level consistency change, not an algorithmic or
precision change on the supported compiler configuration.

## Citation

Please cite the original package paper when using this translation:

Jouni Helske (2017). "KFAS: Exponential Family State Space Models in R."
*Journal of Statistical Software*, 78(10), 1-39.
https://doi.org/10.18637/jss.v078.i10
