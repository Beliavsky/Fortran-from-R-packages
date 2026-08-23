# cbinom-fortran

Modern Fortran translation of the computational core of the R package `cbinom` 1.6 by Dan Dalthorp.

The package implements the continuous analogue of the binomial distribution described by Ilienko (2013), with support `[0, size + 1]`.

## API

- `dcbinom(x, size, prob [, log_p])`
- `pcbinom(q, size, prob [, lower_tail, log_p])`
- `qcbinom(p, size, prob [, lower_tail, log_p])`
- `call rcbinom(x, size, prob)`

Scalar and rank-1 array overloads are provided for d/p/q. `rcbinom` fills a caller-provided rank-1 array.

## Numerical implementation

The upstream package obtains its CDF from R's `pbeta`, numerically differentiates that CDF for the density, and numerically inverts it for quantiles. This port follows the same approach but is self-contained:

- regularized incomplete beta: continued fraction (Lentz-style evaluation)
- density: the same `h = 1e-6` centered/one-sided finite differences used upstream
- quantile: safeguarded monotone bisection to approximately square-root machine precision
- RNG: inverse transform using Fortran `random_number`

For integer `size = n`, `pcbinom(k+1,n,p)` agrees with the ordinary binomial CDF `P(Bin(n,p) <= k)`.

## Build and test

```sh
fpm test
```

or directly with gfortran:

```sh
gfortran -std=f2018 -Wall -Wextra -fcheck=all src/cbinom.f90 test/test_cbinom.f90 -o test_cbinom
./test_cbinom
```

## License

GPL-2.0-or-later, matching the upstream package's `GPL (>= 2)` declaration.
