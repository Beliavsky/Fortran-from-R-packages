# Validation

Validated with GNU Fortran 14.2.0 and the system BLAS/LAPACK libraries.

## Debug / runtime-check build

Library and tests were compiled with:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

All five tests pass:

* `test_compatibility`
* `test_initial_population`
* `test_local_searches`
* `test_malschains`
* `test_rastrigin`

## Optimized build

The same tests and both examples pass with `-O2` and the warning-as-error flags.

Representative release results:

```text
test_malschains:
  fitness = 1.6396e-12
  EA evaluations = 2500
  LS evaluations = 2500

test_rastrigin:
  fitness = 1.6582e-09

example/sphere:
  fitness = 2.8165e-14
```

The 30-dimensional Rastrigin example uses a shorter budget than the upstream R
demo and is included as an executable demonstration rather than a strict
regression threshold.

The linker may emit the standard GNU warning that an executable stack is
required for test/example programs containing internal procedure callbacks.
That comes from the compiler trampoline for the internal Fortran objective
procedure, not from the numerical library.
