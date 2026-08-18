# PoissonBinomial-fortran

Modern Fortran 2018/FPM translation of the computational code in R package
`PoissonBinomial` 1.2.8. It implements ordinary and generalized
Poisson-binomial probability mass functions, CDFs, quantiles and random
number generation without R, Rcpp or FFTW.

## Implemented ordinary methods

Exact methods:

- `DivideFFT` - divide-and-conquer polynomial convolution using a native
  radix-2 FFT for sufficiently large subproblems
- `Convolve` - direct dynamic convolution
- `Characteristic` - characteristic-function DFT inversion
- `Recursive` - the upstream two-column recursive dynamic program

Approximations:

- `Mean`
- `GeoMean`
- `GeoMeanCounter`
- `Poisson`
- `Normal`
- `RefinedNormal`

The ordinary high-level API is:

```fortran
use poisson_binomial, only : dp, dpbinom, dpbinom_at, dpbinom_values
use poisson_binomial, only : ppbinom, ppbinom_at, ppbinom_values
use poisson_binomial, only : qpbinom, qpbinom_values, rpbinom
```

`dpbinom()` and `ppbinom()` return the complete support in order. For a
returned array `a`, `a(k+1)` corresponds to outcome `k`. Scalar and vector
query routines are also provided.

## Generalized Poisson-binomial distribution

The generalized law allows Bernoulli trial `i` to contribute integer
`val_p(i)` with probability `probs(i)` or `val_q(i)` otherwise. Implemented
methods are:

- `DivideFFT`
- `Convolve`
- `Characteristic`
- `Normal`
- `RefinedNormal`

The generalized API uses `type(gpb_table)` for complete tables:

```fortran
type(gpb_table) :: tab

tab = dgpbinom(probs, val_p, val_q, "Convolve")
! Outcome x is tab%values(x - tab%lower + 1).
```

Scalar/vector density and CDF queries, quantiles and RNGs are provided as
`dgpbinom_at`, `dgpbinom_values`, `pgpbinom_at`, `pgpbinom_values`,
`qgpbinom`, `qgpbinom_values` and `rgpbinom`.

Both ordinary and generalized routines accept optional nonnegative integer
multiplicity weights. RNGs support `generator="Sample"` and
`generator="Bernoulli"`.

## Build and test

With FPM:

```text
fpm build
fpm test
fpm run --example example_poisson_binomial
```

The release was independently compiled with GNU Fortran 14.2 using:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all
```

All included tests pass under those flags.

## Files

- `src/` - Fortran library
- `test/` - numerical and RNG tests
- `example/` - small usage example
- `TRANSLATION_NOTES.md` - detailed coverage and implementation differences
- `UPSTREAM.md` - upstream provenance
- `LICENSE` - GNU GPL version 3
- `upstream/` - retained upstream metadata
