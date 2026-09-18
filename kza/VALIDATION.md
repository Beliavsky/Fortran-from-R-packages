# Validation

Validation was performed with GNU Fortran 14.2.0 using Fortran 2018 mode, runtime checking, warnings, implicit-interface warnings, and line truncation promoted to an error.

The maintained modules were compiled directly, in dependency order, with the equivalent of:

```text
gfortran -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=line-truncation -fcheck=all ...
```

The deterministic test executable completed successfully and printed:

```text
All kza tests passed.
```

The test suite exercises every mapped exported computational function, including:

- KZ 1D/2D/3D averaging and missing-value behavior.
- KZA constant-field stability, 4.2.0 full-window semantics, single-slice 3D behavior, matrix rotational symmetrization, max normalization, robust quantile normalization, sparse quantile fallback, and 1D tail marking.
- KZSV on a constant series.
- KZFT zero-frequency values and KZS equivalence.
- KZFT transfer function at its center frequency.
- Raw periodogram transform convention.
- KZ third-order periodogram output shape and finite values.
- Rolling local variance in clamp, zero-pad, `krnl=1`, 2D, and 3D modes.

The example also compiled and ran successfully, printing:

```text
KZ center value:           0.34378
KZA center value:          1.07936
RLV center variance:       1.07301
First spectral magnitude:  30.07769
```

An additional independent reference check compared `kzft` for an even-window case against a direct implementation of the upstream R algorithm; complex values agreed to floating-point roundoff, including the upstream asymmetric padding/slicing convention.

## FPM availability

The requested FPM commands were explicitly attempted from the package root:

```text
fpm build
fpm test
fpm run --example basic_kza
fpm clean --all
```

This sandbox does not contain an `fpm` executable, so each command returned `fpm: command not found`. Direct GNU Fortran compilation, tests, and example execution above are therefore the available build validation and are not represented as an FPM run.

A final static audit checks the maintained Fortran sources and package tree for the repository requirements: free-form `.f90` source only, one package `dp` kind, no `double precision`/`real*8`/D-exponent literals, no semicolon-separated statements, no self-comparison NaN tests, explicit `INTENT` or `VALUE` on every dummy, one dummy per declaration, trailing meaningful FORD `!!` comments on dummy declarations, no duplicate Fortran source files, no vendored dependency source, and internally consistent translation coverage metadata.
