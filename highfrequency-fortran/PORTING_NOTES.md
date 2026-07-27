# Porting notes

## Representation

The R package accepts `xts` and `data.table` objects, dispatches by day and
symbol, and uses named columns. The Fortran library instead accepts explicit
numeric arrays:

- integer time vectors;
- real price, return, size, bid, and offer vectors;
- return matrices with observations in rows and assets in columns;
- typed result structures for jump tests, HAR/HEAVY models, lead-lag results,
  and liquidity measures.

This separates numerical algorithms from storage and presentation.

## Return convention

Upstream `makeReturns` returns an array with the same length as the input and
sets the first return to zero. Fortran `make_returns` returns the mathematically
natural `n-1` differences. Callers that need aligned observations can prepend a
zero explicitly.

## Missing data

No sentinel value is used for missing data. Cleanup routines return logical
masks, and matching routines optionally return availability flags. This avoids
silently treating a NaN as data inside matrix operations.

## Numerical dependencies

The upstream package uses RcppArmadillo, robustbase, Rsolnp, numDeriv,
sandwich, RcppRoll, and other packages. The Fortran translation is
self-contained and includes:

- Gaussian distribution functions and quantiles;
- Gaussian elimination and matrix inversion;
- symmetric Jacobi eigendecomposition and PSD projection;
- ordinary least squares;
- bounded Nelder-Mead optimization.

## HEAVY optimization

The two HEAVY equations are estimated separately with positivity and
stationarity bounds. The upstream package uses Rsolnp and computes robust
sandwich standard errors with numerical derivatives. The Fortran model
preserves the recursions and quasi-likelihood but does not report those standard
errors.

GNU Fortran may emit an executable-stack linker warning because the objective
callbacks are internal procedures. This does not affect the numerical result.
Compilers that implement procedure trampolines differently may not emit it.

## Realized quarticity convention

The Fortran `realized_quarticity` uses the documented formula

```text
n / 3 * sum(r_i^4)
```

for `n` supplied returns. Upstream R data workflows often retain a leading zero
return, which changes the apparent observation count. The Fortran API requires
the caller to supply the intended return sample explicitly.

## Synchronization

- `previous_tick` performs last-observation-carried-forward sampling.
- `refresh_time_pair` implements pairwise refresh-time synchronization.
- `hayashi_yoshida_covariance` works directly on asynchronous intervals.
- Full multivariate refresh-time list construction remains outside the current
  compiled scope.

## Corrections and safeguards

- All covariance outputs are explicitly symmetrized.
- Optional PSD projection floors negative eigenvalues at zero.
- Kernel and regression routines validate dimensions and small-sample cases.
- Fortran logical expressions do not assume short-circuit evaluation.
- Price-based logarithmic routines reject or neutralize nonpositive prices.
- Anticipated zero denominators are handled explicitly rather than propagating
  infinities.
- Cleanup routines return masks instead of mutating caller data.

## Source style

- Modern free-form Fortran 2018
- lower-case source
- `dp = kind(1.0d0)`
- `implicit none`
- standard 132-column line limit
- no external preprocessor
- ASCII translated source
