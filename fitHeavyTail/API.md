# API reference

All public names are available from:

```fortran
use fitheavytail
```

The floating-point kind is `dp = kind(1.0d0)`.

## Result type

The four fitting procedures return a `heavy_tail_fit` derived type through an
`intent(out)` argument. Depending on the model, it contains:

- `mu(:)`: location vector
- `gamma(:)`: skewness vector for `fit_mvst`
- `scatter(:,:)`: scatter estimate
- `covariance(:,:)`: covariance estimate when mathematically defined
- `mean(:)`: distribution mean when mathematically defined
- `nu`: degrees of freedom
- `loadings(:,:)`, `psi(:)`: factor-model covariance decomposition
- `latent_weights(:)`: final conditional scale expectations
- `converged`, `num_iterations`, `cpu_time`, `log_likelihood`
- `status`, `message`

A nonconverged fit may still contain a usable final iterate. Its status is
`ht_no_convergence` rather than `ht_success`.

## `fit_tyler`

```fortran
call fit_tyler(x,result,estimate_mu,initial_mu,initial_covariance, &
   max_iter,ptol,ftol,nu_min,nu_max)
```

Rows containing IEEE NaNs are dropped. The shape matrix is trace-normalized
during iteration, then the package's scale and implied Student-t `nu` recovery
are applied.

## `fit_cauchy`

```fortran
call fit_cauchy(x,result,initial_mu,initial_covariance, &
   max_iter,ptol,ftol,nu_min,nu_max)
```

Implements the accelerated majorization-minimization update from the R source.
Rows containing IEEE NaNs are dropped.

## `fit_mvt`

```fortran
call fit_mvt(x,result,fixed_nu,nu_method,nu_iterative_method, &
   initial_nu,initial_mu,initial_scatter,initial_covariance, &
   na_rm,optimize_mu,observation_weights,scale_covmat, &
   px_em_acceleration,nu_update_start_at_iter, &
   nu_update_every_num_iter,factors,max_iter,ptol,ftol, &
   nu_min,nu_max)
```

Important choices:

- Supply `fixed_nu` to keep `nu` fixed.
- Otherwise `nu_method` defaults to `"iterative"`.
- One-shot methods are `"kurtosis"`, `"cross-cumulants"`,
  `"all-cumulants"`, `"Hill"`, `"MLE-diag"`, and
  `"MLE-diag-resampled"`.
- Iterative methods include `"POP"`, `"POP-approx-1"` through
  `"POP-approx-4"`, `"POP-exact"`, the two sigma-corrected POP methods,
  `"OPP"`, `"OPP-harmonic"`, `"ECM"`, `"ECM-diag"`, `"ECME"`, and
  `"ECME-diag"`.
- `factors < size(x,2)` activates the low-rank plus diagonal factor model.
- `na_rm=.true.` drops incomplete rows. `na_rm=.false.` performs the original
  conditional-moment EM calculations using IEEE NaNs as missing values.
- Observation weights are normalized to have mean one.

## `fit_mvst`

```fortran
call fit_mvst(x,result,fixed_nu,fixed_gamma,initial_nu, &
   initial_mu,initial_gamma,initial_scatter,max_iter,ptol,ftol, &
   pxem,nu_min,nu_max)
```

The generalized-hyperbolic skewed-t estimator drops incomplete rows, as does
the original R function. Omitting `fixed_nu` estimates `nu`; omitting
`fixed_gamma` estimates the skewness vector. The covariance is defined only for
`nu > 4`; otherwise the returned covariance entries are IEEE quiet NaNs.

## Degrees-of-freedom and tail estimators

- `nu_opp_estimator`
- `nu_pop_estimator`
- `nu_mle`
- `excess_kurtosis_unbiased`
- `nu_from_average_marginal_kurtosis`
- `nu_from_cross_cumulants`
- `nu_from_all_cumulants`
- `nu_hill_estimator`
- `nu_pareto_tail_index`

`nu_mle` accepts the internal R method names:

- `MLE-mv-cov`, `MLE-mv-scat`
- `MLE-mv-diagcov`, `MLE-mv-diagcov-resampled`
- `MLE-mv-diagscat`, `MLE-mv-diagscat-resampled`
- `MLE-uv-var-ave`, `MLE-uv-scat-ave`, `MLE-uv-var-stacked`

## Densities and special functions

- `mvt_log_likelihood`
- `mvst_log_likelihood`
- `digamma_dp`
- `log_bessel_k`
- `bessel_k_ratio`
- `sample_skewness`

The Bessel implementation supports real order and uses a scaled integral for
fractional base orders followed by stable log-domain recurrence.

## Status constants

- `ht_success`
- `ht_invalid_argument`
- `ht_too_few_observations`
- `ht_singular_matrix`
- `ht_no_convergence`
- `ht_numerical_error`
