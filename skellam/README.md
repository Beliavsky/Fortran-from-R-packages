# skellam-fortran

A self-contained modern Fortran translation of the computational core of the R
package `skellam` 0.2.4.

The Skellam distribution is the distribution of the difference of two
independent Poisson random variables. This package provides exact and
saddlepoint distribution functions, random generation, maximum-likelihood
estimation, and log-link regression.

## Features

- PMF and log-PMF
- CDF, upper-tail probabilities, and log probabilities
- Quantiles
- Random generation as a difference of Poisson variates
- Saddlepoint PMF and Lugannani-Rice CDF approximations
- Distribution moments
- Two-rate maximum-likelihood estimation
- Skellam regression with separate log links for both Poisson rates
- Numerical covariance matrices, standard errors, Wald statistics, and p-values
- No external numerical-library dependency

## Build with FPM

```text
fpm build
fpm test
fpm run --example distribution_table
fpm run demo_skellam
```

The package version in `fpm.toml` is the FPM-compatible numeric version
`0.2.4`.

## Minimal example

```fortran
program example
   use skellam, only : dp, i8, dskellam, pskellam, qskellam
   implicit none

   print *, dskellam(1_i8, 3.0_dp, 2.0_dp)
   print *, pskellam(1.0_dp, 3.0_dp, 2.0_dp)
   print *, qskellam(0.95_dp, 3.0_dp, 2.0_dp)
end program example
```

## Regression convention

The upstream `skellam.reg()` example models a response formed as
`count2 - count1`, while returning the first coefficient vector for `count1`
and the second for `count2`. `fit_skellam_regression` preserves that convention
by default.

Set

```fortran
response_is_lambda1_minus_lambda2=.true.
```

to use the direct distribution convention `response = count1 - count2`.

Fortran callers pass an `n x p` predictor matrix. An intercept is added by
default and can be disabled with `include_intercept=.false.`.

## Validation

The project contains six standalone regression programs covering distribution
values, quantile inversion, random moments, saddlepoint approximations, MLE,
and regression. See `docs/VALIDATION.md`.

## Scope

Plotting, R vector recycling, R warnings/options, formula processing,
data frames, row/column names, and R list/class infrastructure are not ported.
See `docs/API_MAP.md` and `docs/PORTING_NOTES.md`.

## License

GPL-2.0-or-later. See `LICENSE`, `COPYING`, and `COPYING.GPL-3`.
