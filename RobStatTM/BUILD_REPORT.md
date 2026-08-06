# Build report

Validation date: 2026-08-03

Compiler and libraries:

- GNU Fortran 14.2.0
- System LAPACK and BLAS shared libraries
- Linux x86-64 validation environment
- FPM manifest included; the `fpm` executable was not installed in the validation environment, so the equivalent explicit dependency-ordered builds in `scripts/` were used.

## Checked build

Command:

```text
./scripts/build_checked.sh
```

Key flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface
-fcheck=all -fbacktrace -fimplicit-none -O0 -g
```

Result: all five tests passed.

## Optimized build

Command:

```text
./scripts/build_optimized.sh
```

Key flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface
-Wno-error=maybe-uninitialized -Wno-error=uninitialized
-fimplicit-none -O3 -march=native
```

Result: all five tests passed. GNU Fortran emitted optimizer-only uninitialized-descriptor warnings in the vendored robustbase/rrcov translations. These diagnostics were not observed in the checked build and did not produce runtime failures in either test profile.

## Example output

```text
MM coefficients:      1.00262    2.51833
DCML coefficients:    1.00290    2.50889
DCML LS mixing:       0.01251
```

The synthetic data-generating coefficients were `(1.0, 2.5)` with two large response outliers.
