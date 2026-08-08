# Validation

The release is validated with GNU Fortran using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

The strict test suite contains five programs:

1. `test_rand` - rand/1/bin minimization;
2. `test_current` - current/1/bin minimization;
3. `test_best` - best/1/bin minimization;
4. `test_trace_purify` - non-divisible thinning and periodic purification;
5. `test_options` - defaults and parameter validation.

Two examples exercise the public FPM API.

Run the compiler-independent package tests with:

```text
fpm test
```

or the strict GNU Fortran checks with:

```text
scripts/test_gfortran.sh
```

On Windows:

```text
scripts\test_gfortran.bat
```
