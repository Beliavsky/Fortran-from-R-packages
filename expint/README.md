# expint-fortran

Modern Fortran translation of the computational core of the R package
`expint` 0.2-1 (Vincent Goulet et al.).

The package provides the real exponential integrals

- `E1(x)`
- `E2(x)`
- `En(x)` for integer `n >= 0`
- `Ei(x) = -E1(-x)`

and the upper incomplete gamma function

`Gamma(a,x) = integral_x^infinity t^(a-1) exp(-t) dt`

for real `a` and `x >= 0`, including negative values of `a`.

## Status

Version 0.1.0 translates the numerical kernels in upstream
`src/expint.c` and `src/gamma_inc.c`.  The Chebyshev coefficients and
piecewise algorithms for the exponential integral are direct translations of
the GSL/SLATEC-derived upstream implementation.  The negative-`a` incomplete
gamma continued fraction and recursion follow the upstream code.  For positive
`a`, this Fortran package supplies a self-contained incomplete-gamma
series/continued-fraction implementation in place of the upstream call to
R's `pgamma`.

No external numerical library is required.

## Public API

```fortran
use expint_fortran

real(dp) :: x, y

x = 1.275_dp
y = expint(x)                 ! E1(x)
y = expint(x, 3)              ! E3(x)
y = expint(x, 3, .true.)      ! exp(x) * E3(x)
y = expint_e1(x)
y = expint_e2(x)
y = expint_en(x, 10)
y = expint_ei(x)
y = gammainc(-1.2_dp, 2.5_dp)
```

The scalar functions are `elemental`, so conformable arrays and scalar
broadcasting work naturally in Fortran.

For callers that specifically want R-style recycling of unequal-length
arguments, the package also provides

```fortran
values = expint_recycle(x_array, order_array)
values = gammainc_recycle(a_array, x_array)
```

## Scaling

For `E1`, `E2`, and `En`, `scale=.true.` returns `exp(x) * En(x)`.

As in upstream R `expint_Ei`, the scaled `Ei` interface is implemented as
`-E1(-x, scale=.true.)`; therefore it returns `exp(-x) * Ei(x)`.

## Build with FPM

```text
fpm build
fpm test
fpm run --example basic
```

FPM was not installed in the translation environment, so the complete source,
tests, and example were instead compiled directly with GNU Fortran using

```text
gfortran -std=f2018 -Werror=implicit-interface -Werror=trampolines -fcheck=all
```

## Validation

The permanent tests cover:

- upstream Abramowitz-Stegun reference values for exponential integrals;
- scaled/unscaled identities;
- `E0`, `E1`, `E2`, and higher-order `En`;
- negative arguments accepted by upstream `E1`/`Ei`;
- upper incomplete gamma for positive, zero, negative, and negative-integer `a`;
- the small-`x`, `a` near `-0.5` recursion case highlighted by upstream issue #2;
- invalid-domain NaN behavior;
- R-style argument recycling helpers.

An additional dense validation grid was compared with high-precision `mpmath`.
On that grid, maximum scaled absolute/relative error was approximately
`2.2e-16` for `E1` and `2.6e-13` for `Gamma(a,x)`.

## Differences from the R package

The R `.External` interface, R attributes, warnings, localization, and R object
handling are not translated.  Invalid scalar domains return IEEE NaN, and
underflow returns zero, but the Fortran routines do not emit R warnings.

The native Fortran API uses integer `order`; therefore there is no implicit
truncation of a real-valued order as in R.  Convert explicitly with `int()` if
that behavior is desired.

## License and provenance

The upstream package declares `GPL (>= 2)`, while the numerical C sources
translated here explicitly contain GSL-derived code distributed under GPL
version 3 or later.  This Fortran translation is therefore distributed under
**GPL-3.0-or-later**.

The original `DESCRIPTION`, `CITATION`, R wrappers, and the two translated C
source files are retained in `upstream/` for provenance.  See `COPYING` and
`docs/TRANSLATION_NOTES.md`.
