# CGNM-fortran

Modern Fortran/FPM translation of the computational core of the R package
**CGNM 0.9.3 (Cluster Gauss-Newton Method)**.

The implementation is standalone and uses array/procedure APIs rather than R
functions, environments, data frames, Shiny objects, or plotting classes.

## Implemented

- Cluster Gauss-Newton method, current algorithm version 3 behavior.
- Legacy algorithm-version-1 distance scaling.
- Random or user-supplied initial clusters.
- Independent adaptive Tikhonov regularization for every cluster member.
- Local weighted linear models and regularized CGNR updates.
- K-means neighborhood partitioning and source-style cluster-size fallback.
- Matrix-valued per-member targets and weights.
- Initial-range step restriction.
- Per-parameter lower/upper/both-sided transformations.
- Multi-objective parameter targets (`MO_weights` / `MO_values`).
- `keepInitialDistribution` parameter freezing in the local linear model.
- Residual and lambda histories.
- Residual-resampling, multinomial-weight, and random-weight bootstrap variants.
- EBE weighting/repeated fits.
- `topIndices`, accepted-index screening, best approximate minimizers, and
  column quantiles.
- Low-level CGNR routines corresponding to the R helpers.

## Build

```text
fpm build
fpm test
```

The public module is `cgnm`.

## Small example

```fortran
use cgnm, only : dp, cgnm_problem, cgnm_options, cgnm_result, &
                 cgnm_init_problem, cgnm_fit

type(cgnm_problem) :: prob
type(cgnm_options) :: opt
type(cgnm_result) :: res
integer :: ierr

call cgnm_init_problem(prob, model, target, lower_range, upper_range, ierr=ierr)
opt%num_minimizers = 100
opt%num_iterations = 25
call cgnm_fit(prob, opt, res)
```

See `example/` and `API.md` for complete programs.

## Important differences from R

The package-owned numerical logic is translated, but external R-runtime pieces
are not duplicated. In particular:

- R's `stats::kmeans` is replaced by standalone Lloyd k-means. The cluster
  partition can therefore differ for the same RNG seed even though the CGNM
  local-model logic is unchanged.
- `MASS::ginv` is replaced by a symmetric eigendecomposition pseudoinverse.
- R formula/string reparameterization expressions are replaced by explicit
  bound transformations in `cgnm_problem`.
- The apparent lower-bound sign inconsistency in the R wrapper's internal
  transformation is implemented using the mathematically intended
  `theta = lower + exp(z)` transformation; see `TRANSLATION_COVERAGE.md`.
- Logging to `.RDATA`, Shiny, ggplot2 visualizations, and profile-likelihood
  workflows that depend on saved R objects are not part of the Fortran core.

## License

The upstream package is MIT licensed. The original source tree is retained at
`original/CGNM-master/`.
