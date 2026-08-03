# API

Use the aggregate module:

```fortran
use intrinsicfrp
```

All real calculations use `dp = kind(1.0d0)`. Observations are matrix rows.
Status values are `status_ok`, `status_invalid`, `status_singular`,
`status_nonconverged`, and `status_numerical`.

## Result types

- `vector_result`: `estimate`, `standard_errors`, `selected_indices`, status.
- `screening_result`: SDF coefficients, errors, t statistics, selected factors.
- `hj_result`: squared distance and confidence bounds.
- `rank_test_result`: estimated rank, statistic, p-value, and test sequences.
- `pca_result`: risk premia and selected number of principal components.
- `oracle_control`: Oracle weighting/tuning and inference controls.
- `oracle_result`: penalized premia, errors, scores, and selected penalty.
- `fgx_result`: new-factor SDF coefficients, errors, and selected controls.

## Public package-equivalent procedures

### `TFRP`

```fortran
call tfrp(returns, factors, result, include_standard_errors, hac_prewhite)
```

Computes `Cov(F,R) Var(R)^(-1) E[R]`.

### `FRP`

```fortran
call frp(returns, factors, result, misspecification_robust, &
  include_standard_errors, hac_prewhite, screening_level)
```

`misspecification_robust=.false.` selects Fama-MacBeth; `.true.` selects the
Kan-Robotti-Shanken estimator. A positive `screening_level` first applies GKR
factor screening.

### `SDFCoefficients`

```fortran
call sdf_coefficients(returns, factors, result, misspecification_robust, &
  include_standard_errors, hac_prewhite, screening_level)
```

Computes Fama-MacBeth or GKR linear-factor SDF coefficients.

### `GKRFactorScreening`

```fortran
call gkr_factor_screening(returns, factors, result, target_level, hac_prewhite)
```

Sequentially removes the factor with the smallest insignificant robust SDF
coefficient t statistic.

### `OracleTFRP`

```fortran
call oracle_tfrp(returns, factors, penalties, result, control)
```

`control%weighting_type` is `c` (correlation), `a` (loadings), or `b` (initial
premia). `control%tuning_type` is `g` (GCV), `c` (k-fold CV), or `r` (rolling).

### `FGXFactorsTest`

```fortran
call fgx_factors_test(gross_returns, control_factors, new_factors, result, n_folds)
```

Implements the three-pass factor-selection/test workflow. `gross_returns` follows
the upstream interface and should contain gross, not excess, returns.

### `HJMisspecificationDistance`

```fortran
call hj_misspecification_distance(returns, factors, result, ci_coverage, &
  hac_prewhite)
```

Returns the squared distance and an asymptotic confidence interval.

### Identification tests

```fortran
call iterative_kleibergen_paap_2006_beta_rank_test(returns, factors, result, &
  target_level)
call chen_fang_2019_beta_rank_test(returns, factors, result, n_bootstrap, &
  target_level_kp2006, seed)
```

The Chen-Fang bootstrap is reproducible when a 64-bit `seed` is supplied.

### `GiglioXiu2021RiskPremia`

```fortran
call giglio_xiu_2021_risk_premia(returns, factors, result, which_n_pca, n_max_pca)
```

A positive `which_n_pca` fixes the PCA count. Otherwise the translated selector
chooses it up to `n_max_pca`.

### `HACcovariance`

```fortran
call hac_covariance(series, covariance, status, prewhite)
call hac_variance(series, variance, status, prewhite)
call hac_standard_errors(series, standard_errors, status, prewhite)
```

The input series should already be centered when that is required by the model.

## Lower-level procedures

The aggregate module also exposes moment-form estimators, covariance helpers,
PCA selectors, adaptive weights, Oracle soft thresholding, and FGX covariance
routines. See the module interfaces in `src/` for exact declarations.
