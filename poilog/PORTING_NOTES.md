# Porting notes

## Upstream-to-Fortran mapping

| Upstream routine | Fortran routine | Location |
|---|---|---|
| `dpoilog` / C `poilog` | `dpoilog`, `dpoilog_vec` | `src/poilog_distribution.f90` |
| `dbipoilog` / C `bipoilog` | `dbipoilog`, `dbipoilog_vec` | `src/poilog_distribution.f90` |
| `rpoilog` | `rpoilog` | `src/poilog_rng.f90` |
| `rbipoilog` | `rbipoilog` | `src/poilog_rng.f90` |
| `poilogMLE` | `poilog_mle_fit` | `src/poilog_mle.f90` |
| `bipoilogMLE` | `bipoilog_mle_fit` | `src/poilog_mle.f90` |
| R `optim(..., method="BFGS")` | internal BFGS | `src/poilog_optimize.f90` |
| R `optim(..., method="Nelder-Mead")` | internal Nelder-Mead | `src/poilog_optimize.f90` |
| R/QUADPACK integration | adaptive Gauss-Kronrod 15 | `src/poilog_quadrature.f90` |
| `rnorm` / `rpois` | Box-Muller / PTRS+Knuth Poisson RNG | `src/poilog_rng.f90` |

## Deliberately omitted

- Plotting examples (`barplot`, `hist`, `legend`).
- R data-frame/list presentation details.
- R-specific `write.table` backup output during bootstrap.
- R dynamic-library registration and `.C` interface machinery.

## Intentional robustness changes

- `dbipoilog` requires both standard deviations to be positive. The upstream R
  validation accidentally allows `sig1 == 0`, although the C integrator divides
  by its variance.
- `rho = +/-1` is evaluated through the exact one-latent-normal limit rather
  than passing a zero conditional variance to the nested univariate integral.
- The two apparent R bugs described in `README.md` are corrected.

## Validation

The test suite includes:

- Univariate PMF normalization.
- Independent high-accuracy reference values for `dpoilog` and `dbipoilog`.
- `rho = 0` factorization into the product of univariate PMFs.
- RNG shape/zero-truncation invariants.
- Convergent univariate and bivariate MLE exercises.
