# Porting notes

## Data and precision

All numerical calculations use double precision with
`dp = kind(1.0d0)`. R vectors and lists are represented by assumed-shape arrays
and typed result structures.

The sample skewness and kurtosis definitions match the package: the sample
standard deviation uses denominator `n - 1`, followed by the arithmetic mean of
the standardized third or fourth powers.

## Native optimization

The original package calls `Rsolnp::solnp`. The Fortran port uses a
self-contained constrained Nelder-Mead method. Invalid parameter vectors receive
a large smooth constraint penalty. Numerical Hessians are calculated at the
best solution.

The optimizer is derivative-free and deterministic. The result exposes the
iteration and evaluation counts and returns the best feasible point when the
iteration limit is reached.

## Standard errors

The upstream code calculates:

```text
sqrt(diag(Hessian) / n)
```

This is not the usual likelihood covariance estimator. The Fortran field
`standard_errors` uses `sqrt(diag(inverse(Hessian)))`. The upstream calculation
is retained in `upstream_standard_errors` for comparison.

## AIC and BIC

The upstream functions minimize a negative log likelihood but calculate AIC and
BIC as if that quantity were the ordinary log likelihood. The Fortran fields
`aic` and `bic` use:

```text
AIC = 2*NLL + 2*k
BIC = 2*NLL + k*log(n)
```

The original signed formulas are retained as `upstream_aic` and `upstream_bic`.

## Forecast corrections

Both upstream forecast routines contain this assignment inside the horizon loop:

```text
mu_fcst <- c(h_fcst, ...)
```

It replaces the mean forecast vector with the variance forecast vector. The
Fortran port correctly appends/assigns the next mean forecast.

The upstream code also resets `T` to the length of arrays with the first
observation removed, then indexes the original data with that shortened value.
The Fortran implementation uses the actual final observation and final state.

For GJRSK, multi-step propagation preserves the package's coefficient-sum
approximation for asymmetric terms.

## Initial GJRSK parameters

The upstream initial intercept formulas omit the final persistence coefficient
in each GJRSK block. The Fortran initial values use all dynamic coefficients:

```text
b0 = (1 - b1 - b2 - b3) * variance
c0 = (1 - c1 - c2 - c3) * skewness
d0 = (1 - d1 - d2 - d3) * kurtosis
```

## Likelihood constant

The original active likelihood omits `log(2*pi)`. The Fortran likelihood and
estimation routines match this by default. The optional `include_constant`
argument adds the normal-density constant.

## Parameter constraints

The constraints follow the actual R implementation. In particular, the
skewness intercept `c0` remains unrestricted, while its dynamic coefficients and
their sum are restricted to `(-1, 1)`.
