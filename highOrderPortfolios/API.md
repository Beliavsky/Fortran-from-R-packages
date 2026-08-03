# API

All real arguments use `dp = kind(1.0d0)` from module `fitheavytail_kinds`,
re-exported by `highorderportfolios`.

## Types

### `sample_moments`

Contains `mu`, `covariance`, centered observations, optional `coskewness` and
`cokurtosis` tensors, dimensions, scaling factors, and status information.

### `skew_t_parameters`

Contains skew-t location `mu`, skewness vector `gamma`, `scatter`, its upper
Cholesky factor, degrees of freedom `nu`, and coefficients `a11` through `a43`.

### `portfolio_result`

Contains weights `w`, objective history, four portfolio moments, convergence
status, and, for tilting, `delta` and the four relative improvements.

## Estimation

```fortran
call estimate_sample_moments(x, result, adjust_magnitude, store_tensors)
```

`x` is observations by assets. Covariance uses the usual `T-1` denominator;
third and fourth central co-moments use denominator `T`, consistent with the
upstream portfolio co-moment calculations.

```fortran
call estimate_skew_t(x, result, nu_lb, max_iter, ptol, ftol, pxem)
```

Fits the generalized hyperbolic skew-t model. `nu_lb` must exceed 8 so all four
portfolio moments exist. A `hop_not_converged` result still contains the latest
usable parameter estimates.

## Moment evaluation

```fortran
m = eval_portfolio_moments(w, sample_statistics)
m = eval_portfolio_moments(w, skew_t_statistics)
```

The result is `[mean, variance, skewness, kurtosis]`. Skewness and kurtosis are
third and fourth central moments, not standardized coefficients.

## Sample-moment MVSK design

```fortran
call design_mvsk_portfolio_via_sample_moments(lambda, statistics, result, &
   w_init, leverage, method, tau_w, gamma, zeta, maxiter, ftol, wtol, stopval)
```

Minimizes

```text
-lambda(1) m1 + lambda(2) m2 - lambda(3) m3 + lambda(4) m4
```

on the long-only simplex. Methods: `Q-MVSK`, `MM`, and `DC`. As upstream,
`leverage` currently must equal 1.

## Skew-t MVSK design

```fortran
call design_mvsk_portfolio_via_skew_t(lambda, parameters, result, &
   w_init, method, gamma, zeta, tau_w, beta, tau, initial_eta, &
   maxiter, ftol, wtol, stopval)
```

Methods: `L-MVSK`, `DC`, `Q-MVSK`, `SQUAREM`, `RFPA`, and `PGD`.

## Portfolio tilting

```fortran
call design_mvsktilting_portfolio_via_sample_moments(d, statistics, result, &
   w_init, w0, w0_moments, leverage, kappa, method, gamma, zeta, &
   maxiter, ftol, wtol, theta, stopval)
```

Maximizes the minimum signed, normalized improvement in the four moments while
keeping the tracking variance no larger than `kappa**2`. Methods:
`Q-MVSKT` and `L-MVSKT`.

## Status constants

- `hop_success`
- `hop_invalid_argument`
- `hop_dimension_mismatch`
- `hop_numerical_error`
- `hop_not_converged`
- `hop_fit_error`
