# sandwich-fortran

Modern Fortran 2018 translation of the computational core of the R package
`sandwich` 3.1-2.

The library implements model-robust covariance estimation from estimating
functions (scores) and a bread matrix. It includes heteroscedasticity-consistent,
heteroscedasticity-and-autocorrelation-consistent, clustered, panel, bootstrap,
and jackknife estimators. A self-contained weighted least-squares implementation
is included so the estimators can be used directly with linear models.

## Build with FPM

```text
fpm build
fpm test
fpm run
```

The public module is:

```fortran
use sandwich
```

No external numerical library is required.

## Minimal example

```fortran
use sandwich, only : dp, ols_model, fit_ols, vcov_hc, SANDWICH_SUCCESS
real(dp) :: x(6, 2), y(6)
real(dp), allocatable :: covariance(:, :)
type(ols_model) :: model
integer :: status

x(:, 1) = 1.0_dp
x(:, 2) = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
y = [1.2_dp, 1.9_dp, 3.1_dp, 3.8_dp, 5.2_dp, 5.9_dp]

call fit_ols(x, y, model, status)
if (status /= SANDWICH_SUCCESS) error stop
call vcov_hc(x, model%residuals, model%bread, 'HC3', covariance, status, model%hat)
print *, model%coefficients
print *, covariance
```

## Implemented scope

- Generic bread/meat sandwich covariance and outer-product-of-gradients covariance.
- Weighted OLS fitting, score extraction, bread matrices, and leverage values.
- Constant, HC0, HC1, HC2, HC3, HC4, HC4m, and HC5 covariance estimators.
- Truncated, Bartlett, Parzen, Tukey-Hanning, and quadratic-spectral kernels.
- HAC covariance with explicit weights and VAR prewhitening.
- Andrews AR(1) and ARMA(1,1) bandwidth approximations.
- Newey-West automatic bandwidth and lag-weight generation.
- Lumley-Heagerty WEAVE weights, PAVA, isotonic ACF, and long-run variance.
- One-way and multiway clustered covariance with HC0-HC3 adjustments.
- Panel longitudinal/Driscoll-Kraay and panel-corrected covariance estimators.
- Cluster bootstrap, wild bootstrap, residual bootstrap, fractional bootstrap,
  and jackknife covariance for OLS.
- General covariance and symmetric-matrix numerical utilities required by the
  translated algorithms.

## Interface differences from R

R's `sandwich` package obtains scores and bread matrices through S3 methods for
many model classes. Fortran has no corresponding runtime object system, so this
translation uses explicit arrays. For any model, supply an `n x k` score matrix
and a `k x k` bread matrix. The included `ols_model` type covers weighted linear
least squares directly.

Formula parsing, model frames, `zoo` objects, missing-value actions, R class
attributes, external model-specific adapters, and parallel R callbacks are not
translated. There is no plotting code in the computational package. Bootstrap
refitting is supplied for OLS; other models can use `bootstrap_covariance` or
`jackknife_covariance` with externally generated coefficient replicates.

See `TRANSLATION_COVERAGE.md`, `PORTING.md`, and `API.md` for details.

## License

The upstream package is licensed under GPL version 2 or GPL version 3. This
translation preserves that choice and is distributed under
`GPL-2.0-only OR GPL-3.0-only`. Full texts are in `LICENSE-GPL-2` and
`LICENSE-GPL-3`. The retained upstream snapshot is under `original/`.
