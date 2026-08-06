# StReg modern Fortran

A self-contained modern Fortran translation of the computational code in the
R package **StReg 1.1**, which estimates Student's t static and dynamic
regression models from a joint elliptical distribution.

## Implemented models

- `stlm`: static Student-t linear regression of `y` on contemporaneous `x`
- `star`: Student-t autoregression
- `stdlm`: dynamic Student-t linear regression with contemporaneous exogenous
  variables and lagged endogenous/exogenous variables
- `stvar`: multivariate Student-t vector autoregression

All models retain the upstream block-Toeplitz covariance-factor
parameterization. The fitted joint Student-t distribution is converted to a
conditional mean, a quadratic conditional variance factor, residuals, fitted
values, equation-level diagnostics, and Anderson-Darling tests.

## Build

```sh
fpm test
fpm run --example streg_demo
```

or, without FPM:

```sh
make test
make MODE=optimized test
make example
```

GNU Fortran 14.2 was used for validation. The code requires Fortran 2018 and
no external BLAS, LAPACK, R, or numerical library.

## Basic use

```fortran
use streg, only : dp, streg_fit, streg_options, stlm
real(dp) :: y(100), x(100,2)
type(streg_fit) :: fit
type(streg_options) :: options

options%max_iter = 500
options%compute_hessian = .true.
fit = stlm(y, x, v=6.0_dp, options=options)
```

Fortran arrays use observations by rows and variables by columns. A missing
`trend` argument means an intercept is included. Pass
`include_intercept=.false.` to omit it, or pass a full trend/dummy matrix.

## Main result fields

- `beta`: deterministic coefficients followed by conditional predictor
  coefficients
- `innovation_scale`: conditional Student-t scale matrix
- `conditional_factor`: observation-specific quadratic scale multiplier
- `fitted`, `residuals`, and exact deterministic `trend`
- `joint_scale` and the optimized source parameter vector `theta`
- `r_squared`, `f_statistic`, and `ad_test`
- optional Hessian, transformed-parameter standard errors, and p-values

The fitted conditional covariance for observation `t` is

```text
v / (v + p - 2) * innovation_scale * conditional_factor(t)
```

where `p` is the number of conditioning variables.

## Optimization

`streg_options%optimizer` accepts `"bfgs"` (default) or `"nelder-mead"`.
The upstream sample-size rule is enforced by default; set
`enforce_sample_rule=.false.` only for deliberate small-sample experiments.
The default start is deterministic: trend coefficients are initialized by
least squares and the final block factor is initialized from a symmetric
square root of the residual covariance. A source-format `init` vector can also
be supplied.

## Source compatibility

`source_compatible=.true.` reproduces the upstream residual diagnostic's use
of the response dimension in the conditional Student-t degrees of freedom.
The default uses the mathematically appropriate predictor dimension. Other
source limitations and corrections are detailed in `docs/PORTING_NOTES.md`.

## License

The upstream package declares `GPL-2`. This port is distributed under
GPL-2.0-only to preserve that declaration. See `LICENSE` and `upstream/`.
