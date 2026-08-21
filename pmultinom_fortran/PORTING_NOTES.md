# Porting notes

## Scope

Translated computational routines from `R/pmultinom.R`:

- multinomial probability between lower/upper componentwise bounds
- multinomial CDF and complementary/tail special cases
- inverse one-sided sample-size calculation
- truncated-Poisson convolution machinery, reformulated algebraically as scaled
  unnormalized restricted Poisson convolution

No plotting code exists in the upstream computational file; R documentation,
roxygen, and package/testthat infrastructure were not translated.

## Deliberate implementation changes

1. The upstream package uses the R package `fftw`. This port contains a radix-2
   Cooley-Tukey FFT, so the FPM package has no external numerical dependency.
2. FFT arithmetic uses a real kind with at least 18 decimal digits when the
   compiler provides one, falling back to `dp` otherwise. Public inputs and
   outputs remain double precision (`dp = kind(1.0d0)`).
3. The upstream algorithm divides each restricted Poisson PMF by its truncation
   probability and later multiplies those truncation probabilities back into
   the result. This port cancels those factors symbolically and convolves scaled
   unnormalized restricted Poisson masses. The resulting multinomial
   probability is mathematically identical and avoids fragile tail-probability
   subtraction.
4. R's vector recycling behavior for `lower`, `upper`, and `probs` is retained:
   each shorter vector length must divide the longest vector length exactly.
5. R `NA` behavior has no direct Fortran analogue and is therefore not exposed
   in the Fortran API. Invalid arguments use `error stop`.

## Verification performed

Compiled with GNU Fortran using `-std=f2018 -Wall -Wextra -Wpedantic
-fcheck=all` and ran all test programs successfully. FPM itself was not
installed in the translation environment, so `fpm test` could not be invoked
there; the FPM manifest and directory layout are included for normal FPM use.
