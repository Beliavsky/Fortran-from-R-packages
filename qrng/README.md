# qrng-fortran

Modern free-form Fortran/FPM translation of the computational code in the R
package **qrng 0.0-11** by Marius Hofert and Christiane Lemieux.

The port concentrates on numerical quasi-Monte Carlo functionality and omits
R registration, package-discovery logic, plotting/demo presentation, and
wrappers whose actual algorithms live in other R packages.

## Implemented computational API

The umbrella module is `qrng`.

### Quasi-random generators

- `korobov(n, d, generator, randomize)`
  - scalar or length-`d` generator
  - scalar generators are expanded as in the R interface,
    `generator**(0:d-1) mod n`
  - optional Cranley-Patterson-style random shift used by upstream qrng
- `ghalton(n, d, method, shift_coeff)`
  - generalized Halton and plain Halton
  - all 360 upstream prime bases
  - all Faure-Lemieux generalized-Halton permutation factors
  - the 32-digit randomized shift used by upstream qrng
- `sobol(n, d, randomize, skip, seed, digital_shift)`
  - dimensions through 16,510
  - complete upstream primitive-polynomial and direction-number data
  - skipping of initial points
  - optional 52-bit digital shift

`shift_coeff` and `digital_shift` are optional deterministic controls added to
make parity tests reproducible.  When they are absent, random shifts use
`r_mod`'s RNG helpers.

### Utilities

- `to_array_matrix(x, f)` corresponds to qrng's `format="(n*f,d)"`
- `to_array_3d(x, f)` corresponds to qrng's `format="(n,f,d)"`
- `sum_of_squares(u)`
- `sobol_g(u, alpha)` implements the independent-copula numerical kernel
- `exceedance_indicator(x, q, p)`
- `exceedance_individual_given_sum(x, q, p)`
- `exceedance_sum_given_sum(x, q, p)`

The three exceedance routines are separate because their result ranks/types
differ in R and a type-stable Fortran API is preferable.

## Build

With FPM and BLAS/LAPACK available:

```text
fpm test
fpm run --example example_qrng
```

The full supplied `r_mod` module is retained, so BLAS/LAPACK are linked even
though the qrng routines themselves do not require dense linear algebra.

A direct GNU Fortran validation build used Fortran 2018, preprocessing for
`r_mod.F90`, `-Werror=implicit-interface`, runtime checking, BLAS, and LAPACK.

## Example

```fortran
program example_qrng
use qrng, only: dp, sobol
implicit none
real(dp), allocatable :: u(:,:)
integer :: i

u = sobol(8, 3)
do i = 1, size(u,1)
   print '(3f10.6)', u(i,:)
end do
end program example_qrng
```

## Upstream attribution

The translated generators retain the algorithm provenance stated by qrng:

- Korobov implementation: Marius Hofert, based on C. Lemieux's RandQMC.
- Generalized Halton implementation: Marius Hofert, based on C. Lemieux's
  RandQMC; scrambling factors from Faure and Lemieux (2009).
- Sobol implementation: Marius Hofert, based on Christiane Lemieux's RandQMC.

References cited by the upstream package include:

- Faure, H. and Lemieux, C. (2009), *Generalized Halton Sequences in 2008:
  A Comparative Study*, ACM TOMACS 19(4), Article 15.
- L'Ecuyer, P. and Lemieux, C. (2000), *Variance Reduction via Lattice Rules*.
- Lemieux, C., Cieslak, M. and Luttmer, K. (2004), *RandQMC User's guide*.
- Radovic, Sobol and Tichy (1996), Faure and Lemieux (2009), Owen (2003), and
  Sobol and Asotsky (2003) for the test functions.

See `upstream/DESCRIPTION`, `upstream/*.c`, and `PORTING_NOTES.md` for details.

## License

Code derived from qrng retains the upstream `GPL-2 | GPL-3` license choice.
The supplied `r_mod` module is MIT-licensed and remains separately identified.
