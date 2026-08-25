# pmultinom-fortran

**Official CRAN title:** One-Sided Multinomial Probabilities

Modern Fortran/FPM translation of the computational code in the R package
`pmultinom` 1.0.0 by Alexander Davis. Plotting/document-generation machinery is
not part of the translation.

The package evaluates

`P(lower(i) < N(i) <= upper(i), i=1,...,k)`

for a multinomial count vector, and implements the one-sided inverse sample-size
calculation from the R package.

## Build and test

```text
fpm build
fpm test
fpm run --example wolf_die
```

## Public API

```fortran
use pmultinom_module, only : dp, pmultinom, pmultinom_many, &
                             invert_pmultinom, invert_pmultinom_many
```

Scalar probability:

```fortran
p = pmultinom(size, probs, lower=lower, upper=upper)
```

Either bound may be omitted. Omitted `lower` means `-Inf`; omitted `upper` means
`+Inf`. Bound/probability arrays retain the upstream R package's recycling rule:
shorter lengths must divide the longest length exactly.

Inverse sample size:

```fortran
n = invert_pmultinom(probs, target_prob, lower=lower)
n = invert_pmultinom(probs, target_prob, upper=upper)
```

As in the upstream package, inversion is supported only for one-sided events.

## Numerical implementation

The R package follows Levin (1981): independent Poisson counts with means
proportional to the multinomial cell probabilities are restricted to the desired
bounds, convolved, and then conditioned on their total count.

This port uses an in-tree radix-2 FFT rather than the R package's `fftw`
dependency. It also convolves scaled *unnormalized* restricted Poisson masses.
The upstream code first divides by each Poisson truncation probability and later
multiplies those same factors back. Cancelling those factors algebraically avoids
an unnecessary source of underflow while computing the same conditional
multinomial probability.

## License and provenance

The upstream package declares `License: AGPL-3`. This translation is distributed
under `AGPL-3.0-only`; see `LICENSE`. Original package metadata and the original
computational R source are retained under `upstream/` for provenance and license
preservation.
