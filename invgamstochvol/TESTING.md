# Testing

## FPM

```text
fpm test
```

## Direct GNU Fortran validation

Unix-like systems:

```text
sh scripts/test_gfortran.sh
```

Windows with GNU Fortran available on `PATH`:

```text
scripts\test_gfortran.bat
```

The validation suite builds with:

```text
-std=f2018 -Wall -Wextra -Wpedantic -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow
```

and separately with optimized flags:

```text
-O3 -std=f2018 -Wall -Wextra -Wpedantic -Werror
```

## Test coverage

- hypergeometric identities and a terminating series;
- exact C++ likelihood regression values;
- recursion table dimensions and positivity;
- deterministic exact posterior draws;
- zero persistence and minimal outer truncation;
- invalid lengths, parameters, and nonfinite inputs.
