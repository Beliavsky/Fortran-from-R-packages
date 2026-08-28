# Porting notes

## Upstream basis

This port was made from the attached source tree for **TruncatedNormal 2.3**.
The computational R/C++ sources and package metadata are retained under
`upstream/`.  The translation targets modern standard-conforming free-form
Fortran and the Fortran Package Manager (FPM).

## Computational mapping

| Upstream code | Fortran implementation |
| --- | --- |
| `lnNpr`, `Phinv` | `truncated_normal_math` |
| `norminvp`, `normq`, `qfun` | `truncated_normal_math` |
| `trandn`, `ntail`, `tn`, `trnd` | `truncated_normal_math` |
| `qtnorm`, `rtnorm` | `qtnorm_vec`, `rtnorm_vec` |
| `.cholpermGGE`, `.cholpermGB`, `cholperm` | `truncated_normal_linear:cholperm` |
| `.dmvnorm_arma`, `.dmvt_arma` | `dmvnorm`, `dmvt` |
| `gradpsi`, `jacpsi`, `psy` | `truncated_normal_tilting` |
| `gradpsiT`, `psyT` | `truncated_normal_tilting` |
| `mvnpr`, `mvnprqmc` | same names in `truncated_normal_core` |
| `mvtpr`, `mvtprqmc` | same names in `truncated_normal_core` |
| `mvNcdf`, `mvNqmc` | `mvncdf`, `mvnqmc` |
| `mvTcdf`, `mvTqmc` | `mvtcdf`, `mvtqmc` |
| `mvnrnd`, `mvtrnd` | private proposal kernels in `truncated_normal_core` |
| `mvrandn`, `mvrandt` | same names in `truncated_normal_core` |
| `tregress` | `tregress` returning `tregress_result` |
| `pmvnorm`, `pmvt`, `dtmvnorm`, `dtmvt` | same names in `truncated_normal_api` |
| `ptmvnorm`, `ptmvt`, `rtmvnorm`, `rtmvt` | same names in `truncated_normal_api` |

Fortran is case-insensitive, so names such as `mvNcdf` naturally become
`mvncdf` without changing their callable spelling semantics.

## Reused helpers and dependencies

The port intentionally does not reimplement facilities already provided by the
supplied dependency translations:

- `nleqslv-fortran` solves the minimax saddlepoint systems using the same
  Broyden/Newton trust-region strategy exposed by upstream nleqslv.
- `alabama-fortran` supplies the constrained optimization fallback in the
  Normal tilting solve.
- `spacefillr-fortran` supplies Owen-scrambled Sobol points for `mvnprqmc`.
- `qrng-fortran` supplies randomized digital-shift Sobol points for
  `mvtprqmc`.
- the previously supplied MIT `r_mod.F90`, present in the qrng dependency,
  supplies `normal_cdf`, `qnorm`, Normal/Student RNG/CDF/quantile helpers,
  `solve_real`, and other R-compatible primitives used here.

Only TruncatedNormal-specific missing mathematics was added: robust Normal
interval log probabilities, the Botev truncated-Normal proposal kernels,
Cholesky permutation logic, tilting objectives/derivatives, probability
estimators, and exact samplers.

## Numerical/Fortran adaptations

1. **Observation orientation.** R low-level generators often use a `d x n`
   matrix and transpose in their high-level wrapper.  The Fortran public
   generators consistently return `(n,d)`, matching the high-level R-facing
   convention and simplifying Fortran iteration over samples.

2. **Parameter checks.** R's recycling/default/missing-argument machinery is
   not recreated.  The Fortran APIs use explicit conforming arrays.  Bounds
   used by the low-level minimax routines must have matching dimensions.

3. **Infinite bounds.** IEEE infinities are supported.  `dtmvnorm` and `dtmvt`
   detect completely untruncated regions and skip probability estimation, so
   ordinary multivariate density evaluation is exact and deterministic.

4. **Student normal limits.** The upstream wrapper behavior is retained:
   `dtmvt` and `rtmvt` use the Normal path for `df=0` or `df=+Inf`, and
   `ptmvt` uses the Normal approximation for `df=0` or `df>350`.

5. **Degenerate truncation coordinates.** Upstream high-level generators treat
   widths below `1e-10` as fixed variables and sample the remaining variables
   from the appropriate conditional Normal/Student distribution.  The same
   Schur-complement and Student df/scale update are implemented in the Fortran
   `rtmvnorm` and `rtmvt` wrappers.  `r_mod`'s existing `solve_real` is reused
   for the required solves.

6. **Minimax solver fallback.** The main Normal solve uses `nleqslv`; if it
   does not meet the residual criterion, the port calls the supplied
   `alabama-fortran` augmented-Lagrangian solver rather than implementing a new
   constrained optimizer.  The Student solve follows the upstream Broyden
   path.

7. **QMC dependency split.** This intentionally follows upstream package
   dependencies: the Normal QMC path uses `spacefillr`'s Owen-scrambled Sobol
   generator, while the Student path uses `qrng`'s randomized Sobol generator.

8. **Rcpp/RcppArmadillo.** Their small upstream kernels were mathematical
   density/Cholesky routines, not external algorithms needed at runtime.  They
   are translated to Fortran; no C++ bridge is retained.

## Omitted non-computational material

The port omits R package registration, Rcpp-generated wrappers, S3/UI
presentation, plotting/examples that require R graphics, datasets, vignette
machinery, argument recycling, and package documentation plumbing.  Original
sources are retained only under `upstream/` for provenance and parity review.

## Validation

The current translation was compiled with GNU Fortran in Fortran 2018 mode
with runtime checking and `-Werror=implicit-interface` for the project sources.
The retained tests cover deterministic identities, reference probabilities,
both Cholesky-permutation variants, QMC/MC paths, exact sampling, conditional
fixed-coordinate sampling, and `tregress` reconstruction.

Reference values used for the nontrivial two-dimensional checks:

- Normal rectangle: `0.44175424561021764`
- Student-t (`df=5`) rectangle: `0.39766111423138145`

For `d=15` with unit marginal variances and all off-diagonal correlations
`0.5`, the positive-orthant Normal probability is `1/(d+1)=1/16`; this is also
an upstream tinytest target and is retained here.

FPM was not available in the execution environment used for the translation,
so validation used an equivalent direct GNU Fortran build against the vendored
dependency sources plus BLAS/LAPACK.  The included `fpm.toml` is structured for
normal `fpm build` / `fpm test` use.
