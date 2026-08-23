# greybox-fortran

Modern Fortran/FPM translation of the numerical core of R package `greybox` 2.0.8.

Version 0.2 expands the initial distribution/regression port with the major remaining standalone modeling workflows: beta regression, occurrence/hurdle models, ARIMA-error ALM, robust/penalized losses, dynamic point-IC model averaging, automatic demand identification, RMCB and distribution-shape-robust bootstrap.

## Main capabilities

### Custom distributions

Complete d/p/q/r implementations are provided for Laplace, asymmetric Laplace, generalized normal/exponential-power, S, folded normal, Box-Cox normal, logit-normal, rectified normal and three-parameter lognormal.

### Forecast measures and diagnostics

ME, MAE, MSE, MRE, MIS, MPE, MAPE, MASE, RMSSE, relative/scaled measures, GMRAE, SAME, pinball loss, half moments, asymmetry and extremity are available.

### Matrix/data helpers

The port includes polynomial/backshift utilities, lag/lead expansion, regressor products, transformations, dynamic multipliers, outlier/temporal dummies, Cramer's V, multiple/partial correlation and determination diagnostics.

### ALM regression

`alm_fit` supports the main continuous, count and binary families from upstream greybox, response-scale prediction, point likelihoods, covariance estimates and AIC/AICc/BIC/BICc.

v0.2 additionally provides:

- beta regression with separate shape-1 and shape-2 coefficient blocks;
- LASSO, RIDGE, ROLE and QUALE fitting losses;
- `alm_fit_occurrence` for logistic/probit occurrence plus a positive-demand model;
- `alm_fit_arima_errors` for conditional ARIMA-error `orders=(p,d,q)` fitting.

### Model selection and combination

- `stepwise_fit`
- `calm_fit` / `calm_predict`, including beta-regression model averaging
- `lm_dynamic_fit`: observation-specific point-IC model weights with LOWESS-style smoothing
- `rolling_origin_alm`
- `recursive_lm`
- coefficient bootstrap and scale-model regression

### Demand and forecast-comparison workflows

- `aid_fit` and `aid_cat`
- `rmcb_test` with Tukey/Nemenyi, normal-rank and ALM branches
- `dsr_bootstrap` with additive/multiplicative, parametric/nonparametric and intermittent-demand modes

## Build

With FPM:

```sh
fpm test
fpm run --example example_alm
```

The source layout was validated directly with GNU Fortran using:

```sh
gfortran -std=f2018 -Wall -Wextra -Werror -fcheck=all
```

No external numerical library is required.

## Validation

The accumulated suite contains:

```text
test_alm: PASS
test_core: PASS
test_parity_v02: PASS
test_selection: PASS
test_stats_utils: PASS
```

`test_parity_v02` exercises beta ALM and beta model averaging, all four newly added loss modes, occurrence models, AR-error fitting, Tukey and non-Gaussian RMCB branches, standard/intermittent bootstrap, AID/AIDCat and dynamic point-IC combination.

## Scope

The library is matrix-native Fortran. R formula/model-frame handling, S3 methods, plotting (`spread`, `tableplot`, `graphmaker`), `zoo` indexing and xtable/texreg presentation are intentionally excluded.

See `API_MAP.md` and `PORTING_NOTES.md` for exact mappings and the remaining lower-priority differences.
