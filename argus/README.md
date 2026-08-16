# argus-fortran

Modern Fortran implementation of the computational code in the R package
`argus` 0.1.1.

The library implements the Argus distribution on `0 <= x <= 1`, `chi > 0`:

    f(x) proportional to x * sqrt(1-x^2) * exp(-chi^2*(1-x^2)/2)

It has no R or Runuran dependency.

## Features

- density and log-density: `dargus`
- lower/upper CDF: `pargus`
- log CDF/tail probabilities
- lower/upper quantiles: `qargus`
- log-probability quantile input
- inversion random generation
- ratio-of-uniforms random generation
- varying-parameter random generation
- elemental array evaluation
- explicit R-style recycling helpers
- standalone shape-3/2 incomplete-gamma calculations
- FPM project layout

## Build with FPM

```text
fpm build
fpm test
fpm run --example basic
```

## Basic use

```fortran
program demo
   use argus, only : dp, dargus, pargus, qargus
   implicit none

   print *, dargus(0.3_dp, 1.0_dp)
   print *, pargus(0.3_dp, 1.0_dp)
   print *, qargus(0.9_dp, 2.0_dp)
end program demo
```

Expected values are approximately:

```text
0.728912103796695
0.109495296387555
0.939855232066375
```

## Random generation

```fortran
use argus, only : dp, rargus, seed_argus_rng, ARGUS_INVERSION, ARGUS_ROU
real(dp) :: x(10000)

call seed_argus_rng(12345)
call rargus(x, 0.3_dp, ARGUS_INVERSION)
call rargus(x, 3.0_dp, ARGUS_ROU)
```

For one parameter per observation use `rargus_varying`.

## Vector evaluation

The distribution functions are elemental:

```fortran
real(dp) :: x(3), f(3)
x = [0.1_dp, 0.5_dp, 0.9_dp]
f = dargus(x, 0.3_dp)
```

For R-style recycling of unequal numeric argument lengths, use
`dargus_recycle`, `pargus_recycle`, or `qargus_recycle`.

## License

GPL-2.0-or-later, matching the upstream package. See `LICENSE`,
`LICENSES.md`, `NOTICE.md`, and `upstream/`.
