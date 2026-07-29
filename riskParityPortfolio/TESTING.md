# Testing

The port was tested with GNU Fortran 14.2.0 and system LAPACK/BLAS.

## Strict debug build

The library, all tests, the application, and all examples were compiled with:

```text
-std=f2018 -O0 -g -Wall -Wextra -Werror -Wimplicit-interface
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

The following test programs pass:

- `test_core`: published eight-asset weights, equalized risk budgets, agreement
  among Spinu, Roncalli, Choi, and Newton solvers, active-risk limiting case,
  diagonal analytical solution, and high-level API.
- `test_formulations`: finite-difference checks for every objective gradient and
  every Jacobian column across all eight formulations, plus SCA solutions.
- `test_constraints`: exact budget-box projection, KKT residuals, general linear
  projection, group constraints, additional equalities, bounds, mean-return
  reward, and variance penalty.

The published `pyrb` reference portfolio is reproduced to better than `2e-6` in
weights and `1e-9` in relative risk contributions.

## Optimized build

An optimized clean build is also run with `-O3 -DNDEBUG`, while retaining
warnings as errors. All tests and runnable targets must pass before packaging.

## FPM validation

The project uses the standard FPM directory layout. The execution environment
used to create this port did not contain an `fpm` executable, so the manifest is
parsed with Python's TOML parser and targets are compiled directly with
`gfortran`. The requested link libraries are listed in `fpm.toml`.
