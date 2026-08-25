# pbivnorm-fortran

**Official CRAN title:** Vectorized Bivariate Normal CDF

Modern Fortran/FPM translation of the R package `pbivnorm` 0.6.0.

The package evaluates the CDF of a standard bivariate normal distribution using
Alan Genz and Yihong Ge's quadrature algorithm from the upstream package.

## API

```fortran
use pbivnorm_mod, only : dp, pbivnorm
real(dp) :: p
p = pbivnorm(0.5_dp, 1.0_dp, 0.4_dp)
```

`pbivnorm` is elemental. Scalar arguments broadcast naturally over conformable
arrays, for example:

```fortran
p = pbivnorm(x, y, 0.4_dp)
```

For R-style recycling of vectors of arbitrary different lengths, use
`pbivnorm_recycle`.

## Build

```text
fpm build
fpm test
fpm run --example example_pbivnorm
```

No external numerical dependencies are required.

## License

GPL-2.0-or-later, matching upstream `License: GPL (>= 2)`.
