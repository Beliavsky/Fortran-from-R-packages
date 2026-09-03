# MCMCglmm

Modern free-form Fortran translation of the computational code of the R package
**MCMCglmm 2.36** by Jarrod Hadfield.

The project is designed to live as a top-level sibling directory in
`Beliavsky/Fortran-from-R-packages`. It preserves upstream GPL licensing and
provenance while replacing R object plumbing and the bundled modified CSparse
implementation with explicit modern Fortran APIs and shared repository
linear-algebra/tree packages.

## Implemented computational core

The public `mcmcglmm` module includes:

- deterministic state-passed random number generation and MCMC variate samplers;
- pedigree relationship/inverse matrices, pruning, and breeding-value simulation;
- phylogenetic precision and breeding values via the sibling `ape` translation;
- ordinary and conditioned inverse-Wishart sampling;
- truncated/conditional multivariate-normal calculations;
- multivariate Gaussian mixed-model Gibbs MCMC;
- multiple independent random-effect `G` structures;
- heterogeneous Gaussian/non-Gaussian multi-`G` latent MCMC with a full `R`;
- grouped multinomial, `ztmb`, zero-truncated multinomial, and two-process
  zero-inflated/hurdle/zero-altered multi-`G` engines;
- a typed `mcmcglmm_fit_numeric` orchestration entry point for heterogeneous
  scalar, parameter-expanded scalar, two-process, grouped multinomial,
  ordinal, threshold, `theta_scale`, structural path/SIR, and `covu` models,
  accepting either dense or canonical CSR `X`/`Z` inputs and using one
  normalized covariance-routing configuration where applicable;
- a one-based canonical CSR matrix type with coordinate/dense construction,
  validation, conversion, products, transpose products, and crossproducts;
- all computationally active upstream family-code likelihood/response kernels;
- distinct ordered-probit and threshold samplers with adaptive cutpoints;
- native binary slice updates and MCMCglmm-style adaptive proposal recursions;
- missing-response imputation across maintained latent family engines;
- covariance update modes 0-6, including antedependence and direct-sum modes;
- Gaussian and heterogeneous multi-`G` parameter-expanded MCMC;
- `theta_scale` Gaussian multi-term MCMC;
- discrete-Berkson `me()` category updating integrated into the heterogeneous
  scalar orchestration path;
- structural path/SIR matrices and an integrated Jacobian-aware Gaussian sampler;
- specialized Gaussian `covu` joint random/residual covariance MCMC;
- posterior covariance transformations and `posterior.mode` KDE;
- numerical `buildV`, posterior linear prediction, response-scale expectations,
  and posterior-predictive latent/response simulation;
- numerical cores corresponding to `Tri2M`, `Ptensor`, `kunif`, `KPPM`, `knorm`,
  `krzanowski.test`, `gelman.prior`, `Ddivergence`, `mult.memb`, `path`, `sir`,
  and the implemented LRTP `spl` basis.

See `API_COVERAGE.md` for the exact family/update-code map and remaining gaps.

## Dependencies

No translated dependency, CSparse, BLAS, LAPACK, or ARPACK source is vendored in
this package. The FPM manifest uses sibling path dependencies:

```toml
rfortran-core = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
ape = { path = "../ape" }
```

`rfortran-core` supplies the shared `dp` kind. `rfortran-linalg` supplies dense
factorizations/solves/eigenanalysis through the repository's pinned pure-Fortran
LAPACK dependency. The translated `ape` package supplies the shared phylogenetic
tree type and tree algorithms.

There are no `-lblas`, `-llapack`, or system-library links in this package.

## Build

With the sibling dependencies present at the repository root:

```text
cd MCMCglmm
fpm build
fpm test
fpm run --example '*'
```

All maintained source is free-form `.f90` and uses the shared `dp` real kind.

## Examples

The `example/` directory contains deterministic programs covering:

- Gaussian MCMC;
- family/ordered engines;
- threshold MCMC;
- parameter expansion;
- heterogeneous multi-`G` MCMC;
- grouped multi-`G` MCMC;
- `theta_scale`;
- structural/path MCMC;
- `covu` joint G-R MCMC;
- `buildV`/prediction;
- pedigree and phylogenetic calculations;
- distribution and tensor utilities.

Every sampler accepts an explicit `rng_state`, allowing reproducible independent
chains without hidden global RNG state.

## Scope and remaining work

Upstream MCMCglmm's production sampler is a large C++/C sparse latent-variable
engine with extensive R formula/object machinery. This translation now contains
most of the numerical pieces and several integrated dense engines, including
multi-`G`, grouped, special covariance, parameter-expansion, theta-scale,
path/SIR, measurement-error, and `covu` paths.

It does **not** claim a drop-in replacement for R `MCMCglmm()`. The
`mcmcglmm_fit_numeric` API is the first common orchestration layer: callers supply
numeric model matrices through `mcmcglmm_numeric_model`, priors through
`mcmcglmm_numeric_prior`, and MCMC/covariance routing through `mcmcglmm_control`.
Dense and canonical CSR design matrices are accepted. The heterogeneous scalar
chain keeps paired CSR `X`/`Z` designs sparse while forming predictors and normal
equations, assembles the joint coefficient precision in CSR storage, and uses an
independently implemented sparse Cholesky factorization and triangular solves.
The stacked design crossproduct is also accumulated directly as CSR. A
column-assisted row accumulator consolidates repeated contributions immediately,
so its storage is proportional to the sparse inputs, canonical output, and column
workspace rather than all observation-level pair contributions. A deterministic
reverse Cuthill-McKee ordering reduces sensitivity to input column order. Other
engine routes materialize CSR inputs before entering their dense kernels. The
result is tagged by engine and contains the corresponding retained chain.

The largest remaining target is extending this entry point into a genuinely joint
sparse engine that can mix multiple response-block classes and special design
routes in one chain with scaling comparable to upstream CSparse. The scalar CSR
route caches its ordering and symbolic Cholesky analysis, rebuilding safely when
the precision graph expands. It still needs more sophisticated ordering,
production-scale benchmarking, and sparse execution in the other engine routes.
Discrete Berkson measurement-error categories are integrated into the
heterogeneous scalar route, but the current specialized routes retain the
combinations supported by their underlying kernels. R formula parsing,
factors/contrasts, `newdata`, S3/coda objects, and plotting remain interface code
and are intentionally omitted.

## Licensing and provenance

MCMCglmm declares `GPL (>= 2)`. This translation is distributed under
GPL-2.0-or-later. See `COPYING`, `NOTICE.md`, and the byte-preserved upstream
metadata under `upstream/`.

For the upstream scientific citation, see `upstream/CITATION`.
