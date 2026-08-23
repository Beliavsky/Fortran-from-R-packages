# ChernoffDist-fortran

Modern Fortran translation of the computational code in CRAN package
`ChernoffDist` 0.1.0.

The package evaluates Chernoff's distribution for

`Z = argmax_t { B(t) - t^2 }`,

using the Groeneboom-Wellner numerical representation used by the upstream R
package.

## Implemented API

```fortran
use chernoffdist, only: dp, dchern, pchern, qchern
```

- `dchern(x [, log_value])`: density
- `pchern(q [, lower_tail, log_p])`: distribution function
- `qchern(p [, lower_tail, log_p])`: quantile
- `dchern_vector`, `pchern_vector`, `qchern_vector`: array convenience wrappers

The three upstream computational exports (`dChern`, `pChern`, `qChern`) are all
covered.

## Numerical implementation

The density follows the upstream two-representation algorithm exactly in
structure:

- the first 20 small-`t` series coefficients;
- the first 20 Airy zeros in the large-`t` auxiliary series;
- the first 100 Airy zeros and `Ai'(a_k)` values in the negative-argument
  representation;
- adaptive 15-point Gauss-Kronrod quadrature for the improper-integral terms;
- symmetry for the distribution function.

The first 100 Airy zeros and derivatives are embedded as binary64 constants, so
the Fortran library has no runtime dependency on R, the R `gsl` package, or a
system GNU GSL installation.

## Build with FPM

```text
fpm test
fpm run --example basic
```

There are no external package dependencies.

## Example

```fortran
program basic
  use chernoffdist, only: dp, dchern, pchern, qchern
  implicit none

  print '(a,f14.10)', 'density at 0 = ', dchern(0.0_dp)
  print '(a,f14.10)', 'cdf at 1     = ', pchern(1.0_dp)
  print '(a,f14.10)', 'q(.90)       = ', qchern(0.90_dp)
end program basic
```

Expected output is approximately

```text
density at 0 =   0.7583445581
cdf at 1     =   0.9752206566
q(.90)       =   0.6642351962
```

## Validation

The source tree is tested with

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

Tests include independent numerical reference values, symmetry, numerical
verification that the CDF derivative agrees with the density, quantile/CDF
inversion, and option/vector interfaces.

## License

The upstream package declares `GPL-3`; the translated code is distributed under
GPL-3.0-only. See `LICENSE`, `LICENSES.md`, and `UPSTREAM.md`.
