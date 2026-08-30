# classInt-fortran

Modern Fortran translation of the computational functionality of the R package
`classInt` 0.4-11.  The project uses the Fortran Package Manager (FPM) layout
and preserves the upstream GPL-2.0-or-later licensing and provenance.

## Implemented functionality

The public module is:

```fortran
use classint
```

The main entry point is `classint_fit`.  It supports the class interval styles:

- `fixed`
- `sd`
- `equal`
- `pretty`
- `quantile` (R quantile types 1 through 9)
- `kmeans`
- `hclust`
- `bclust`
- `fisher`
- `jenks`
- `dpih`
- `headtails`
- `maximum`
- `box`

Additional native APIs include:

- `classify_intervals` and `find_cols`
- `class_counts`
- `get_hclust_class_intervals`
- `get_bclust_class_intervals`
- `fisher_exact`
- `jenks_breaks`
- `dpih_bandwidth`
- `pretty_breaks`
- `jenks_tests`, `gvf`, `tai`, and `oai`
- `classint_loglik`, `classint_aic`, and `classint_n_partitions`

`classify_intervals` can either classify an existing `class_intervals` object or
fit and classify observations in one call.

## Example

```fortran
program demo
    use classint
    implicit none

    type(class_intervals) :: fit
    real(dp) :: x(10)
    integer, allocatable :: cls(:)

    x = [2.0_dp, 3.0_dp, 3.5_dp, 4.0_dp, 8.0_dp, &
         9.0_dp, 10.0_dp, 17.0_dp, 18.0_dp, 30.0_dp]

    call classint_fit(x, 4, "fisher", fit)
    call classify_intervals(fit, cls)
end program demo
```

See `example/classint_example.f90` for a complete program.

## Dependencies

`bclust` uses the translated `e1071-fortran` implementation.  That dependency,
and its `proxy-fortran` dependency, are vendored as FPM path dependencies so
this release is self-contained.

The direct plug-in bandwidth calculation used by style `dpih` is a native
Fortran translation of the computational algorithm in `KernSmooth::dpih` and
its binned kernel-functional estimator.  Its Hermite/Gaussian kernel is wrapped,
zero-padded to the same power-of-two length, and convolved with the binned counts
using a native radix-2 FFT, matching KernSmooth's computational structure.

## Numerical and parity notes

The original `classInt/src/fish1.f` Fisher dynamic-programming kernel has been
retained under `upstream/` and independently compiled during development.  A
deterministic fixture gives the same optimal partitions, class minima/maxima,
means, and population standard deviations to floating-point roundoff.

The saved upstream `classInt` regression value
`logLik(classIntervals(rep(1:3, each=10), n=2, style="jenks")) = -14.52876`
is also covered by the Fortran tests.

There are a few documented language/runtime adaptations:

- random sampling, k-means starts, and bagged clustering use the deterministic
  native RNG supplied by `e1071-fortran`; streams are not bit-identical to R's
  RNG;
- `hclust` is specialized to one-dimensional data and supports `complete`,
  `single`, `average`, `centroid`, `ward.d`, and `ward.d2`; uncommon R linkage
  variants and exact R tie-ordering are not claimed;
- `pretty` follows R's default 1/2/5/10 unit-selection behavior for ordinary
  finite ranges, but exact behavior at extreme subnormal/overflow-scale ranges
  is not claimed;
- `dpih` now follows KernSmooth's FFT convolution layout; last-bit differences
  can still occur because the native radix-2 transform does not use R's FFT
  implementation or its exact floating-point evaluation order;
- date/time, `units`, factor/shingle, S3 printing, colors, and plotting are R
  interface or presentation features and are intentionally not reproduced.

## Source rules

Maintained Fortran source follows the repository translation rules:

- one shared public `dp = real64` kind;
- `real(dp)` variables and `_dp` real constants;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- every dummy argument has explicit `INTENT` or `VALUE`;
- every dummy argument has its own declaration line;
- every dummy declaration has a meaningful trailing FORD `!!` comment;
- free-form source remains within the standard 132-column limit and is
  `fprettify`-compatible.

`fprettify` was not installed in the translation environment, so it was not
run on this release.

## Building

With FPM installed:

```text
fpm build
fpm test
fpm run --example classint_example
```

A compiler-only strict validation harness is also provided:

```text
tools/run_strict_tests.sh
```

It compiles the vendored dependencies and this package using bounds checking,
implicit-interface diagnostics, and floating-point traps before running the
tests and source-rule audits.

## License

The translated package is distributed under GPL-2.0-or-later, following the
upstream `classInt` package.  See `NOTICE.md`, `licenses/`, and `upstream/`.
