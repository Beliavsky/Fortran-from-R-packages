# goftest-fortran

Modern Fortran/FPM translation of the computational core of the R package
`goftest` 1.2-3.

The library implements the null distributions and tests for the classical
Anderson-Darling and Cramer-von Mises goodness-of-fit statistics. It is
self-contained and does not require R or an external special-functions library.

## Main API

```fortran
use goftest, only : dp, p_ad, q_ad, p_cvm, q_cvm
use goftest, only : ad_test, cvm_test, ad_test_values, cvm_test_values
use goftest, only : gof_result, recognise_cdf
```

`p_ad`, `q_ad`, `p_cvm`, and `q_cvm` accept either scalars or rank-1 arrays.
Omitting `n` selects the asymptotic null distribution. For Anderson-Darling,
`fast=.false.` selects the slower high-accuracy Marsaglia asymptotic algorithm.

For data already transformed through the null CDF, use `ad_test_values` or
`cvm_test_values`. For arbitrary data, `ad_test` and `cvm_test` accept a scalar
CDF callback:

```fortran
result = ad_test(x, my_cdf)
```

Set `estimated=.true.` to request Braun's adjustment for a composite null.
The routine generates the same nearly-balanced random grouping structure as the
R implementation. An optional integer `groups(:)` argument can be supplied for
reproducible grouping.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The source is standard modern Fortran and has also been validated directly with
GNU Fortran using strict warnings and run-time checking.

## Numerical implementation

- Anderson-Darling finite-sample and asymptotic calculations are ports of the
  Marsaglia algorithms shipped in the R package.
- The high-accuracy asymptotic Anderson-Darling calculation uses the standard
  `erfc` intrinsic in the place explicitly allowed by the original C comments.
- Cramer-von Mises probabilities use the Csorgo-Faraway expansion from the R
  package.
- Fractional modified Bessel functions required by that expansion are evaluated
  internally from the defining integral, with an asymptotic branch for large
  arguments.
- Quantiles use bracketed bisection rather than R's `uniroot`.

See `PORTING_NOTES.md` and `docs/TRANSLATION_COVERAGE.md` for details.

## License

GPL-2.0-or-later, matching the upstream `License: GPL (>= 2)` declaration.
