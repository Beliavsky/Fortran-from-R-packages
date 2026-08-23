# Porting notes

## Source basis

This port was made from the supplied `rangen-master.zip`, version 0.0.1 (published 2026-03-16). The package contains C++17/Rcpp code and header-only numerical implementations.

All 48 generated R entry points have a numerical Fortran counterpart or direct API equivalent. Rcpp object conversion and OpenMP orchestration are not required by the Fortran library.

## RNG architecture

The upstream package uses PCG32 for uniform/integer generation and the external `zigg` package for normal generation. The supplied archive does not contain `zigg`'s Ziggurat source.

The Fortran port therefore:

1. preserves the PCG32 XSH-RR recurrence,
2. evaluates the 64-bit recurrence modulo 2^64 using 16-bit limbs, so no nonstandard signed-overflow assumption is needed,
3. uses independent PCG32 streams for ordinary distribution uniforms, direct `runif`, normals, and sampling,
4. uses Box-Muller for normal generation.

Consequently, distributional behavior is preserved but exact normal-dependent bit streams are not expected to match the R package.

## Corrections to upstream behavior

### `Rcauchy`

The upstream `Cauchy::operator()` returns

```text
location + scale * tan(pi * norm())
```

where `norm()` is a standard-normal draw. The documentation and the Cauchy inverse-CDF formula require a uniform draw. The Fortran port uses

```text
location + scale * tan(pi * (U - 0.5))
```

with `U` uniform on (0,1).

### Sampling with replacement

The upstream integer generator uses

```text
pcg32_random_r() % max_bound + min_bound
```

For bounds `0, n-1`, this uses modulus `n-1`, so the final element cannot be selected. The Fortran implementation uses modulus equal to the actual range width.

### `Sample.int(n, size < n)`

The upstream helper allocates a result of length `size`, then writes the values `1:n` into it before sampling. That is out of bounds when `size < n`. The Fortran implementation constructs a pool of length `n` and performs a partial Fisher-Yates draw.

### Gamma shape below one

The upstream Marsaglia-Tsang implementation reuses the final acceptance uniform as the power-transform uniform for `shape < 1`. The standard boosted-shape algorithm requires a separate uniform draw. The Fortran implementation uses the independent transform draw, and tests explicitly cover subunit gamma and beta shapes.

### Seeding

In upstream PCG objects, `setSeed()` changes `state` but not the stream increment `inc`, which was initialized from wall-clock time. Thus the same explicit seed need not identify the same PCG stream across process starts. The Fortran API makes seeding deterministic by using fixed, distinct odd stream increments.

## Uniform endpoints

The public `runif` interface retains the upstream closed `[min,max]` conversion from a 32-bit integer. Internal inverse-transform generators use an open `(0,1)` conversion to prevent accidental `log(0)` or endpoint infinities.

## Column sampling

The R package allows a different requested sample length by column/row, but its Armadillo assignment requires matching target dimensions in practice. The Fortran routines allocate to the maximum requested length and fill unused cells with IEEE NaN.

The `parallel` and `cores` optional arguments are accepted by `col_sample` and `row_sample` for source-level familiarity, but v0.1 executes these routines serially. This affects performance only, not the sampling definition.

## Timer

R's `nanoTime()` returns system-clock nanoseconds since the platform epoch through C++ `system_clock`. Standard Fortran does not define an epoch for `system_clock`; `nano_time()` therefore returns nanosecond-scaled `system_clock` ticks suitable for elapsed-time measurement, not a portable Unix-epoch timestamp.
