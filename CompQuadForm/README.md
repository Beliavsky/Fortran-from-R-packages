# CompQuadForm-fortran

Modern Fortran/FPM translation of the computational core of the R package
CompQuadForm 1.4.4.

## Implemented routines

- `davies` - Davies AS 155 survival probability, diagnostic trace, and `ifault`
- `farebrother` - Farebrother/Ruben AS 204 survival probability, density, and `ifault`
- `imhof` - Imhof survival probability and numerical integration error estimate
- `liu` - Liu-Tang-Zhang moment-matching approximation

The public routines use the same core arguments as the R package: `q`,
`lambda`, integer degrees/multiplicities `h`, and noncentralities `delta`, with
matching defaults. Davies also exposes `sigma`, `lim`, and `acc`; Farebrother
exposes `maxit`, `eps`, and `mode`; Imhof exposes `epsabs`, `epsrel`, and
`limit`.

## Build

```text
fpm build
fpm test
fpm run --example demo
```

BLAS and LAPACK are linked because the supplied general-purpose `r_mod.f90`
contains procedures that reference them, even though CompQuadForm's own four
algorithms do not require BLAS/LAPACK directly.

## Example

```fortran
use r_mod, only: dp
use compquadform_mod, only: davies_result_t, davies

type(davies_result_t) :: fit
real(dp) :: lambda(3), delta(3)
integer :: h(3)

lambda = [0.5_dp, 1.2_dp, 2.0_dp]
h = [1, 2, 3]
delta = [0.0_dp, 0.7_dp, 0.3_dp]
fit = davies(6.0_dp, lambda, h, delta)
print *, fit%qq, fit%ifault
```

`qq` corresponds to CompQuadForm's `Qq`, i.e. the upper-tail/survival
probability.

## Licensing and attribution

Translated CompQuadForm-derived code is GPL-2.0-or-later. The supplied
`r_mod.f90` is MIT-licensed. See `LICENSE`, `LICENSES/`, `NOTICE`, and
`upstream/`.
