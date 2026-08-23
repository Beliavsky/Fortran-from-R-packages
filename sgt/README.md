# sgt-fortran

Modern Fortran/FPM translation of the computational code in Carter Davis's
`sgt` R package 2.0.

The library implements the skewed generalized t distribution and the limiting
skewed generalized-error and uniform cases, together with maximum-likelihood
fitting.

## Implemented API

Use the umbrella module:

```fortran
use sgt
```

The public numerical API includes:

- `dsgt`
- `psgt`
- `qsgt`
- `rsgt`
- `sgt_pdf`, `sgt_logpdf`, `sgt_cdf`, `sgt_quantile`
- `sgt_mle_constant`
- `sgt_mle_model`
- `sgt_params`
- `sgt_mle_result`

`dsgt`, `psgt`, and `qsgt` are elemental, so scalar parameters can be applied
to arrays directly.

## Distribution conventions

The defaults reproduce a standard normal distribution:

```fortran
real(dp) :: f
real(dp) :: inf
inf = ieee_value(0.0_dp, ieee_positive_inf)
f = dsgt(1.2_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, inf)
```

As in the R package:

- `q = Inf` gives the skewed generalized-error distribution.
- `p = Inf` gives a uniform distribution.
- `mean_cent=.true.` makes `mu` the mean when it exists.
- `var_adj=.true.` rescales the internal scale so the resulting standard
  deviation is `sigma` when the variance exists.
- `sigma_multiplier=k` together with `var_adj=.false.` corresponds to the R
  package's numeric `var.adj=k` behavior.

Invalid parameter combinations return IEEE NaN.

## Maximum likelihood

`sgt_mle_constant` is the direct replacement for fitting one constant set of
SGT parameters to a sample. Any subset of `mu`, `sigma`, `lambda`, `p`, and `q`
can be held fixed.

`sgt_mle_model` is the array/callback replacement for the R package's formula
interface. A user callback maps an optimization vector and observation index
to the transformed observation and its `mu`, `sigma`, `lambda`, `p`, and `q`.
This supports regression specifications without embedding an R formula parser
in the numerical library.

The default fit runs both Nelder-Mead and BFGS and retains the fit with the
higher likelihood, matching the default strategy of `sgt.mle`. Numerical
scores, Hessians, covariance matrices, standard errors, z statistics, and
normal-approximation p-values are returned in `sgt_mle_result`.

## Build

```text
fpm build
fpm test
```

The source has also been validated directly with GNU Fortran using:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

## Tests

The test suite covers independent reference values, finite-SGT/SGED/uniform
branches, normal/Cauchy/Laplace identities, CDF-quantile inversion, tail and
log options, RNG moments, constant-parameter MLE, and callback regression MLE.

## Licensing

The translated `sgt` code is GPL-3.0-or-later, matching upstream
`License: GPL (>= 3)`.

The supplied `optimx-fortran` and `numDeriv-fortran` archives are included only
under `reference-dependencies/`. They are not linked into the FPM library; see
`LICENSES.md` and `PORTING_NOTES.md` for the conservative license rationale.
