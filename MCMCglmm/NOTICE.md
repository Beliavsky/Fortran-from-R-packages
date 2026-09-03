# NOTICE and provenance

## Upstream package

This project is a modern Fortran translation of computational code from:

- **R package:** MCMCglmm
- **Upstream version:** 2.36
- **Upstream date:** 2024-05-06
- **Author and maintainer:** Jarrod Hadfield
- **Upstream repository:** `https://github.com/jarrodhadfield/MCMCglmm`
- **Upstream license:** GPL (>= 2)

The uploaded source archive used for this translation was
`MCMCglmm-master.zip`. Byte-preserved upstream `DESCRIPTION`, `NAMESPACE`, and
`inst/CITATION` files are included under `upstream/`.

The primary upstream citation is Jarrod D. Hadfield (2010), "MCMC Methods for
Multi-Response Generalized Linear Mixed Models: The MCMCglmm R Package",
*Journal of Statistical Software* 33(2), 1-22. The exact upstream citation file
is preserved in `upstream/CITATION`.

## Translated source families

The maintained Fortran source is newly written, but computational behavior and
algorithms are derived from upstream MCMCglmm code, including in particular:

- `src/MCMCglmm.cc`: mixed-model coefficient and covariance conditionals;
  family-code likelihood/latent updates; missing-response handling; multiple
  random-effect blocks; redundant-parameter (`alpha`) updates; native binary
  slice branches; adaptive acceptance recursions; `theta_scale`; structural
  Lambda calculations; and `covu` joint random/residual updates;
- `src/dcutpoints.c`: ordered/threshold cutpoint likelihood and proposal-ratio
  calculations used by adaptive cutpoint samplers;
- `src/inverseA.cc` and its R wrapper: pedigree relationship/inverse operations;
- `src/rbv.cc`: pedigree breeding-value recursion;
- `src/rtnorm.c`: truncated-normal simulation;
- `src/rtcmvnorm.c`, `src/pcmvnorm.c`, and `src/cs_dcmvnorm.c`: conditional and
  truncated multivariate-normal calculations;
- `src/rIW.cc`, `src/cs_rinvwishart.c`, and `src/cs_rCinvwishart.c`:
  inverse-Wishart and conditioned inverse-Wishart algorithms;
- `src/cs_rR.c`, `src/cs_rRsubinvwishart.c`, and `src/cs_rSinvwishart.c`:
  correlation-constrained, submatrix-constrained, and direct-sum covariance
  updates;
- `src/rante.cc` and `src/cs_rAnte.c`: antedependence covariance calculations;
- `src/cs_schur.c`: Schur-complement decomposition used by `covu`;
- `src/pkk.c`: `pkk` probability calculation;
- `R/MCMCglmm.R`: family/update-code setup, special-term routing, parameter
  expansion, theta-scale, measurement-error, path/SIR, `mev`, and `covu`
  bookkeeping used to interpret the native numerical kernels;
- `R/me.R`: discrete Berkson measurement-error prior-category construction;
- `R/path.R` and `R/sir.R`: path/SIR design-matrix calculations;
- `R/buildV.R`: marginal covariance assembly for prediction;
- `R/predict.MCMCglmm.R`: linear-predictor marginalization and response-scale
  posterior means;
- `R/simulate.MCMCglmm.R`: posterior-predictive random-effect/residual draws and
  family response transformations;
- `R/posterior.mode.R`: posterior-mode KDE behavior;
- R-level numerical routines corresponding to `posterior.cor`,
  `posterior.inverse`, `posterior.evals`, `posterior.ante`, `Tri2M`, `Ptensor`,
  `kunif`, `KPPM`, `knorm`, `krzanowski.test`, `gelman.prior`, `Ddivergence`,
  `prunePed`, `mult.memb`, and `spl`.

The phylogenetic branches of upstream `inverseA` and `rbv` depend on `ape`.
This translation deliberately reuses the sibling translated `ape` package
rather than copying tree algorithms into MCMCglmm.

## CSparse attribution

Upstream MCMCglmm includes a modified CSparse codebase. Its `src/cs.h` identifies
CSparse version 2.2.1 and states:

- Copyright (c) Timothy A. Davis, 2006-2007.

No CSparse source, header, object, or library has been copied into this
translation. This package has independently implemented canonical CSR storage,
natural-order sparse Cholesky factorization, and triangular solves. Dense routes
continue to use shared factorization from the sibling `rfortran-linalg` package.
This notice retains provenance for algorithms that were originally expressed
through upstream sparse operations without vendoring that dependency.

## Intentional interface substitutions

- R's global RNG is replaced by an explicit deterministic `rng_state` passed to
  every stochastic maintained API.
- R formulas/model matrices are replaced by already-built numeric design arrays.
- Modified CSparse storage/solves are replaced by independently implemented CSR
  design operations and sparse Cholesky, plus dense factorization from sibling
  `rfortran-linalg` for routes not yet migrated.
- The specialized `covu` sampler follows upstream's numerical restriction that
  a `ginverse` is not combined with `covu`; prediction/simulation for `covu` is
  also not invented because upstream methods reject that combination.

## Translation license

Because the translated work is derived from GPL-licensed upstream code, the
maintained translation is licensed under GPL-2.0-or-later. The full GPL version
2 license text is provided in `COPYING`; the "or later" grant follows the
upstream `GPL (>= 2)` declaration.
