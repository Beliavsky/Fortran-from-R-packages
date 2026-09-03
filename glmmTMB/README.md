# glmmTMB — modern Fortran computational translation

This directory contains a modern free-form Fortran translation of the portable
computational kernels in the R package **glmmTMB 1.1.14**.  The upstream package
fits generalized linear mixed models with TMB, including zero inflation,
multiple dispersion models, many response families, and structured Gaussian
random effects.

The translation focuses on numerical functionality that can be exposed as a
standalone Fortran API.  R formula processing, S3/S4 methods, plotting,
printing, data-frame manipulation, R serialization, OpenMP/R runtime glue, and
other R-specific interfaces are intentionally excluded.

## Build

Place this directory at repository top level next to the translated `TMB`
package:

```text
Fortran-from-R-packages/
  TMB/
  glmmTMB/
```

Then from `glmmTMB/` run:

```text
fpm build
fpm test
fpm run --example basic_glmmtmb
```

The FPM dependency is a sibling path dependency.  No TMB source is copied into
this package and no system BLAS/LAPACK library is linked.

## Public API

Use the umbrella module:

```fortran
use glmmtmb
```

Major translated functionality includes:

- glmmTMB family, link, covariance-structure, and prior integer codes;
- stable forward/inverse link transformations and zero-truncation probabilities;
- beta-binomial, generalized-Poisson, skew-normal, Bell, mean-parameterized
  COM-Poisson, and compound-Poisson Tweedie density kernels;
- Wishart, inverse-Wishart, and LKJ correlation-prior kernels from the upstream
  C++ distribution helpers;
- conditional observation log likelihoods for the families implemented in
  `src/glmmTMB.cpp`, including zero inflation and zero truncation;
- variance functions corresponding to the computational family definitions;
- Gaussian random-effect negative log likelihoods for diagonal, unstructured,
  compound-symmetry, AR(1), heterogeneous AR(1), Ornstein-Uhlenbeck,
  exponential, Gaussian, Matérn, Toeplitz, reduced-rank, proportional, and
  equal-to covariance structures;
- dense linear-predictor construction and a portable joint negative
  log-likelihood composition API;
- normal, gamma, Student-t, Cauchy, and LKJ prior kernels.

See `API_COVERAGE.md` for exact scope and known differences.

## Numerical notes

The package uses a single package-local `dp = real64` kind.  The sibling TMB
translation also uses `real64`, so its explicit interfaces are kind-compatible.
No maintained source uses `double precision`, `real*8`, D-exponent literals,
or self-comparison NaN tests.

The Matérn correlation implementation evaluates the integral representation of
`K_nu(x)` numerically because standard Fortran has no intrinsic modified Bessel
K function of arbitrary real order.  The deterministic test suite checks it
against an independent reference value.

The COM-Poisson implementation solves for lambda from the exact requested mean
and evaluates its normalizing constant by a finite log-scale sum.  Very large
means that would require more than roughly 100,000 support points return NaN
rather than silently producing a poor approximation.

## Full glmmTMB fitting versus the portable kernel

Upstream glmmTMB performs automatic differentiation through TMB and integrates
Gaussian random effects with a Laplace approximation while optimizing the
marginal likelihood.  That AD/Laplace optimizer runtime is not reproduced here.
The translated `glmmtmb_joint_nll` is the joint objective before Laplace
integration and is suitable for callers that provide their own optimization,
differentiation, or integration strategy.

This distinction is important: the package translates the reusable numerical
model kernels, but it is not a drop-in replacement for the R `glmmTMB()` fitting
front end.

## Validation

`test/test_glmmtmb.f90` contains deterministic checks against independently
computed reference values.  See `BUILD_VALIDATION.md` for commands actually
run in the translation environment and its FPM/fprettify limitation.

## License and provenance

The upstream package declares **AGPL-3**.  This translation is distributed
under **GNU AGPL version 3 only**.  `LICENSE` contains the license text.
`upstream/DESCRIPTION`, `upstream/CITATION`, and `upstream/MD5` are preserved
from the supplied upstream snapshot.  See `NOTICE.md` for provenance details.
