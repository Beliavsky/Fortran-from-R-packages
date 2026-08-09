# Validation

The translation was compiled with GNU Fortran using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

## Regression programs

1. `test_reference_values.f90`
   - evaluates all 50 functions at `lower + 0.37*(upper-lower)`
   - compares against values obtained by compiling and calling the original
     `src/objFun.c`
   - validates the dispatcher and all bound metadata simultaneously
2. `test_known_optima.f90`
   - checks exact/easy known optima for Ackley, Becker-Lago, Goldstein-Price,
     Griewank, Rastrigin, Rosenbrock, Salomon and Schaffer1
3. `test_metadata.f90`
   - checks all 50 dimensions, bounds and published global optimum metadata
   - checks unknown-name and wrong-dimension error paths

All three regression programs pass.

Both examples also compile and run successfully under the same flags.

FPM was not installed in the validation environment, so the FPM source tree was
compiled directly with gfortran.  The manifest is plain TOML and the source
layout follows the normal `src/`, `test/`, and `example/` FPM conventions.
