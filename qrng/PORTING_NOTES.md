# Porting notes

## Source basis

This port was made from the attached `qrng-master.zip`, whose DESCRIPTION
identifies qrng version 0.0-11, authors Marius Hofert and Christiane Lemieux,
and license `GPL-2 | GPL-3`.

The attached `r_mod.f90` was supplied as MIT-licensed and is reused rather
than duplicating applicable R-compatible helper procedures.

## What was translated

### Native C algorithms

The following native qrng implementations were translated directly:

1. `src/korobov.c`
   - Korobov lattice sequence
   - scalar-generator expansion is performed in the Fortran interface
   - optional random shift

2. `src/ghalton.c`
   - plain and generalized Halton sequence
   - all 360 prime bases
   - all 360 Faure-Lemieux permutation factors
   - 32 base-digit random shift

3. `src/sobol.c`
   - Sobol sequence
   - `skip` handling through Gray-code state
   - 52-bit digital shift
   - complete 16,510 primitive-polynomial table
   - complete 16,509 x 17 initial direction-number table

The Sobol and Halton numerical tables were mechanically extracted from the
upstream C source.  A full-table parity check compared every translated entry:
16,510 primitive polynomials, 280,653 initial direction numbers, 360 primes,
and 360 generalized-Halton factors.

### R-level computational helpers

Translated:

- `to_array()` as rank-specific `to_array_matrix()` and `to_array_3d()`
- `sum_of_squares()`
- the numerical independent-copula kernel of `sobol_g()`
- all three numerical modes of `exceedance()` as separate type-stable routines

For the empirical thresholds in `exceedance()`, qrng explicitly requests
`quantile(..., type=1)`.  The supplied `r_mod` `quantile()` currently always
implements Type 7 even when a `type` argument is supplied.  Therefore this
port adds one genuinely missing local helper, `empirical_quantile_type1()`,
instead of silently changing qrng's threshold definition.

## Deliberately not reimplemented

The following are not algorithms implemented by qrng itself and were not
copied from external packages:

- Sobol `randomize="Owen"`, which delegates to
  `spacefillr::generate_sobol_owen_set()`.
- Sobol `randomize="Faure.Tezuka"` and `"Owen.Faure.Tezuka"`, which delegate
  to `randtoolbox::sobol()`.
- The general-copula transform inside `sobol_g()`, which delegates to
  `copula::cCopula()`.

For a non-independent copula, callers can apply the corresponding inverse
Rosenblatt transform externally and pass the transformed uniforms to the
Fortran `sobol_g()` kernel.

R `.Call` registration (`src/init.c`), R package namespace/loading code,
examples whose purpose is plotting/presentation, and package-installation
checks were also omitted.

## API differences made for Fortran

- R's single `to_array()` returns either rank 2 or rank 3 according to a
  string.  Fortran uses two type/rank-stable functions.
- R's `exceedance()` can return a logical vector, matrix, or numeric vector.
  Fortran exposes three separate functions.
- `ghalton()` accepts optional `shift_coeff(d,32)` for deterministic testing.
- `sobol()` accepts optional `digital_shift(d)` for deterministic testing and
  an optional integer `seed` convenience argument.
- `korobov()` is overloaded for scalar and vector generators.

These changes preserve the numerical algorithms while avoiding dynamic
R-style return types.

## r_mod handling

`upstream/r_mod-original.f90` is the attached file verbatim.  The build copy
is `src/r_mod.F90`; only formatting/encoding changes were made:

- long free-form lines were wrapped with continuation markers so no compiler
  option for unlimited line length is needed;
- the uppercase `.F90` suffix lets common Fortran compilers preprocess the
  existing conditional compilation directives automatically;
- the UTF-8 BOM was removed.

No qrng helper already supplied by `r_mod` was independently reimplemented.

## Numerical validation

The retained test programs exercise:

- Korobov scalar-generator expansion and point values;
- plain Halton values with a fixed zero digital shift;
- generalized Halton values and Faure-Lemieux permutation factors;
- Sobol reference points in dimensions 1-4;
- Sobol `skip` behavior;
- Sobol dimension 16,510;
- deterministic randomized-shift reproducibility;
- array reshaping semantics;
- sum-of-squares and Sobol-g kernels;
- Type-1 empirical exceedance thresholds and all exceedance return forms.

A first randomized Sobol build exposed an important Fortran porting issue:
`floor()` without a `kind=` argument returns a default integer, which overflowed
when converting the upstream 52-bit digital shift.  The final code uses
`floor(..., kind=int64)`, matching the upstream 64-bit integer intent.

All retained tests pass under GNU Fortran 14.2 with Fortran 2018,
`-Werror=implicit-interface`, and runtime checking.  FPM itself was not
installed in the execution environment, so the manifest was validated but
`fpm test` could not literally be invoked here.
