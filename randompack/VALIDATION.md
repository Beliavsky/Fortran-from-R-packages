# Validation

This file records validation of the current translation checkpoint.

## Performed in this environment

Compiler: GNU Fortran 14.2.0.

The maintained Fortran sources were compiled with:

```text
-std=f2018 -Wall -Wextra -Werror=line-truncation -fcheck=all
```

Because the sibling `rfortran-linalg` checkout is not present in this sandbox,
a temporary, external test-only `r_linalg` stub exposing the single interface used
by `randompack` (`cholesky_factor`) was used for this direct compiler check. The stub is not included in this package and is not a
substitute for the required integration build against the real sibling package.

Results:

- all maintained source files compiled successfully;
- `test/test_randompack.f90` ran successfully and printed
  `All randompack tests passed.`;
- `example/basic_randompack.f90` compiled and ran successfully;
- the example produced five finite normal draws;
- no line exceeded the standard 132-column free-form limit;
- maintained Fortran source contains no forbidden `double precision`, `real*8`,
  `kind(0.0d0)`, or D-exponent real literals;
- semicolons occur only in the explicitly permitted simple `SELECT CASE` form;
- no self-comparison NaN test is used;
- every dummy argument has explicit `INTENT` or `VALUE`, is declared separately,
  and has a trailing FORD `!!` documentation comment;
- no duplicate maintained Fortran source file was found;
- no dependency source has been copied into the maintained implementation;
- the PCG64 seed-123 discrete sequence was cross-checked directly against the
  upstream C implementation: 64 `int(-5,17)` draws, a 25-element permutation,
  a 10-of-40 sample, and the following 31 raw bytes match exactly after the
  documented one-based index conversion for `perm` and `sample`;
- the corrected ranlux++ first four seed-123 words were cross-checked directly
  against upstream C and frozen into the engine regression vector.
- PCG64 seed-123 normal and exponential sampling was cross-checked directly
  against the upstream C implementation: the first ten outputs of each method
  match bit-for-bit and are frozen into the deterministic test suite;
- forced rare Ziggurat tail paths were also cross-checked against upstream C: normal seed `9071` matches bits `400EF9C97BE4D17D`, and exponential seed `753` matches bits `4020DC22F1CF89A2`;
- the translated double-precision BSD/OpenLibm `log`, `log1p`, and `exp` routines match upstream C bit patterns on ordinary, near-boundary, subnormal, and large reference arguments;
- PCG64 seed-123 gamma vectors were cross-checked directly against upstream C: 17 draws at `shape=2.5, scale=1.2` and 17 draws at `shape=0.7, scale=2.0` match bit-for-bit in `bitexact` mode;
- PCG64 seed-123 beta, Student-t, and F vectors were cross-checked directly against upstream C: 17 draws from `Beta(2,5)`, `t(7)`, and `F(5,11)` match bit-for-bit in `bitexact` mode;
- PCG64 seed-123 transformed-family vectors were cross-checked directly against upstream C: 17 draws each from `lognormal(mu=0.3,sigma=1.4)`, `gumbel(mu=0.2,beta=1.3)`, `pareto(xm=1.7,alpha=2.3)`, `weibull(shape=1.7,scale=2.1)`, and `skew_normal(mu=0.4,sigma=1.2,alpha=-2.5)` match bit-for-bit in `bitexact` mode;
- PCG64 seed-123 multivariate-normal generation was cross-checked directly against upstream C for a 5-by-3 SPD covariance and a rank-2 covariance that forces complete pivoting: all 15 values in each output matrix match bit-for-bit, and the six uniforms immediately following each call also match, verifying RNG consumption order;

The exact-real-equality warnings in three tests are intentional: those tests
verify exact stream restoration/duplication behavior rather than approximate
floating-point numerical equality.

## FPM limitation

`fpm` is not installed in this execution environment. Consequently the requested
release commands

```text
fpm build
fpm test
fpm clean --all
```

could not be executed here. An end-to-end FPM build against the actual sibling
`../rfortran-linalg` dependency therefore remains a release validation item.
The package manifest is configured for that sibling path dependency.

This ZIP should be treated as the current working checkpoint rather than a claim
that the unavailable FPM integration test has passed.
