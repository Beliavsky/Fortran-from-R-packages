# Validation

Validated with GNU Fortran 14.2.0 on Linux.

## Bounds-checked build

Compiler options:

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Werror -Wimplicit-interface
```

Passing regression executables:

- `test_discrete_lp`
- `test_functions`
- `test_geometry`
- `test_global_optimizers`
- `test_history`
- `test_local_optimizers`

## Optimized build

The same six tests pass with:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface
```

Representative example output in this environment:

```text
multiple-knapsack value = 345.00
Nelder-Mead Rosenbrock f ~= 7e-13
simpleDE Rastrigin f = 0
```

The multiple-knapsack MILP can have more than one optimal bin assignment; the
regression checks objective and capacity feasibility rather than requiring one
specific optimal assignment.

GNU ld emits an executable-stack warning for binaries involving internal
procedure callbacks.  This arises from GNU Fortran callback trampolines and is
not a compiler or numerical test failure.

FPM was not installed in the validation environment.  Both `fpm.toml` files
were parsed as TOML and the FPM-equivalent source/dependency build was compiled
directly with gfortran.
