# vars

Modern free-form Fortran translation of the computational algorithms in the R package
`vars` 1.6-1 by Bernhard Pfaff, with contributions by Matthieu Stigler.

The library provides VAR estimation and lag selection, restrictions, forecasting,
MA/orthogonalized MA representations, roots, impulse responses, FEVD, residual
diagnostics, Granger and instantaneous causality tests, Blanchard-Quah
identification, SVAR scoring, and array-level VECM/SVEC calculations.

This directory is intended to live at the repository root beside the shared
`rfortran-core` and `rfortran-linalg` packages.

## Build

```text
fpm build
fpm test
```

The FPM manifest uses sibling path dependencies:

```toml
rfortran-core = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

No source from those packages is copied into this package.

## Numerical kind

All maintained Fortran source uses the single `dp` kind from
`rfortran-core`:

```fortran
use r_kinds, only : dp
```

The public `vars` module re-exports `dp`. Real variables use `real(dp)` and
real literals use the `_dp` suffix.

## Main API

The public module is `vars`.

- `fit_var`: estimate an unrestricted VAR by least squares.
- `var_select`: AIC, HQ, SC/BIC, and FPE lag-order selection using the fixed
  sample convention of the R package.
- `restrict_var_manual`, `restrict_var_ser`: manual and sequential-elimination
  coefficient restrictions.
- `var_loglik`: Gaussian reduced-form VAR log likelihood.
- `phi_from_a`, `psi_from_a_sigma`: MA and orthogonalized MA coefficient arrays.
- `var_roots`: companion-matrix roots.
- `forecast_var`, `forecast_covariance`: point forecasts and forecast-error
  covariance matrices.
- `impulse_response`, `fevd_var`: reduced-form IRF and FEVD calculations.
- `structural_impulse_response`, `structural_fevd`: IRF and FEVD from a supplied
  structural impact matrix.
- `bq_identification`: Blanchard-Quah long-run identification.
- `arch_test_univariate`, `arch_test_multivariate`: ARCH-LM diagnostics.
- `jarque_bera_univariate`, `jarque_bera_multivariate`: normality diagnostics.
- `portmanteau_tests`, `bg_serial_tests`: multivariate residual serial-correlation
  diagnostics.
- `granger_causality`, `instantaneous_causality`: causality tests.
- `svar_negloglik`, `svar_fit_scoring`: structural VAR likelihood and scoring
  estimator under A/B zero restrictions.
- `vec2var_coefficients`: convert VECM coefficient arrays to a level-VAR lag
  representation for either upstream `transitory` or `longrun` specification.
- `svec_long_run_matrix`, `svec_fit_scoring`: array-level SVEC long-run multiplier
  and scoring estimator.
- `residual_bootstrap_path`, `bootstrap_irf_indices`: deterministic residual
  bootstrap calculations using caller-supplied resampling indices.

See `API_COVERAGE.md` for a file-by-file parity map and intentional omissions.

## Minimal example

```fortran
program demo
   use r_kinds, only : dp
   use vars, only : var_model, fit_var, var_const, vars_success
   implicit none

   type(var_model) :: model
   real(dp) :: y(8, 2)
   integer :: info

   y(:, 1) = [1.0_dp, 1.2_dp, 1.1_dp, 1.4_dp, 1.6_dp, 1.5_dp, 1.8_dp, 2.0_dp]
   y(:, 2) = [0.4_dp, 0.5_dp, 0.7_dp, 0.6_dp, 0.9_dp, 1.0_dp, 1.1_dp, 1.3_dp]

   call fit_var(y, 1, var_const, model, info)
   if (info /= vars_success) error stop "fit failed"

   write (*, '(a)') "Coefficient matrix:"
   write (*, '(*(f12.6,1x))') model%coef
end program demo
```

Two complete examples are in `example/`.

## Array conventions

For a `K`-variable VAR(`p`), lag matrices are stored as
`a(K,K,p)`. `a(:,:,j)` multiplies the observation at lag `j`. MA/IRF arrays
use the third dimension for horizon and include horizon zero as element 1.
The coefficient matrix in `var_model%coef` has one row per response equation;
lagged endogenous regressors precede deterministic, seasonal, and exogenous
regressors.

## Scope

The translation intentionally excludes R formula handling, `lm`/`mlm`/S3 object
assembly, print/summary/plot/fanchart methods, and model callbacks. The
`stability()` method is not duplicated because upstream `vars` delegates its
actual computation to `strucchange::efp`; use the translated top-level
`strucchange` package for that computation.

For SVEC/VECM work, this package accepts numerical coefficient arrays rather than
an R `ca.jo` object. This avoids copying the translated `urca` package and keeps
the numerical boundary explicit.

## License and provenance

The upstream R package declares `GPL (>= 2)`. This translation is distributed as
`GPL-2.0-or-later`; both GPL-2.0 and GPL-3.0 texts are included. See `NOTICE.md`
and `UPSTREAM.md` for authorship, source-archive identity, attribution, algorithm
mapping, and citations.
