# ADGofTest modern Fortran port

This package translates the computational code of R package `ADGofTest` 0.3
to modern Fortran. It implements the Anderson-Darling goodness-of-fit statistic
and the finite-sample distribution approximation from Marsaglia and Marsaglia
(2004).

## Features

- Anderson-Darling statistic for a sample already transformed to uniform values.
- Goodness-of-fit testing through a user-supplied scalar CDF procedure.
- Marsaglia finite-sample CDF and upper-tail p-value approximation.
- Source-compatible treatment of exact probabilities 0 and 1: infinite
  statistic and p-value 0.
- Optional endpoint clipping for CDF underflow or overflow.
- Optional clamping of approximate CDFs and p-values to `[0, 1]`.
- Explicit result type, status codes, FPM manifest, Makefile, example, and tests.

R-specific `htest` objects, expression capture, and ellipsis argument dispatch are
not reproduced. Distribution parameters can be captured in a module procedure
or wrapper CDF supplied to `ad_test`.

## Build

With FPM:

```text
fpm test
fpm run --example adgoftest_demo
```

With GNU Make:

```text
make check
make release
```

On Windows with GNU Fortran, run `scripts\\test.bat`.

## Basic use

```fortran
use adgoftest, only : dp, ad_test, ad_test_result

type(ad_test_result) :: result
real(dp) :: x(100)

call ad_test(x, standard_normal_cdf, result, clamp_p_value=.true.)
```

The callback must have the interface

```fortran
pure function standard_normal_cdf(x) result(probability)
   real(dp), intent(in) :: x
   real(dp) :: probability
end function standard_normal_cdf
```

## License

The upstream package declares `License: GPL` without selecting a GPL version.
This port preserves that declaration. See `LICENSE` and `licenses/`.
