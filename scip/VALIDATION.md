# Validation

Validation was performed against the exact SCIP 10.0.2 / SoPlex 8.0.2 source
bundled in the uploaded `scip` package.  The backend was compiled locally with
CMake and linked to the translated Fortran API.

## Strict Fortran build

The library and tests compile with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

The C/C++ ABI shim also compiles warning-clean with `-Wall -Wextra -Werror`.

## Regression tests

Five executables pass:

1. `test_one_shot` — dense and CSC LP equivalence plus binary knapsack.
2. `test_model_api` — incremental variables/linear constraints, objective
   sense, native parameters, solution pool, and solve statistics.
3. `test_special_constraints` — quadratic, SOS1, SOS2, and indicator handlers.
4. `test_status_control` — unbounded/infeasible statuses and one-shot controls.
5. `test_sos_weights` — explicit SOS2 ordering weights.

The production LP from the package vignette gives minimized objective `-22`
at `(10/3, 4/3)`, i.e. maximized profit `22`.

Representative special-constraint checks:

```text
quadratic: min x subject to x^2 >= 4, x >= 0  -> x ~= 2
SOS1:      max x1+x2, 0<=x<=1                -> objective 1
SOS2:      max x1+x2+x3                       -> objective 2
indicator: max 100*z+x, z=>x<=2               -> (z,x)=(1,2), objective 102
```

The complete test set also passes with the translated Fortran compiled at
`-O2`.

FPM itself was not installed in the validation container.  `fpm.toml` is a
normal FPM manifest; the convenience build scripts construct the vendored
backend and place it on `LIBRARY_PATH` before invoking FPM.
