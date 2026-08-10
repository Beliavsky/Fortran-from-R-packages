# Validation

The release is validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

## Permanent regression programs

1. `test_algorithms.f90`
   - Misra1a benchmark.
   - LM, LM+geodesic acceleration, dogleg, double dogleg, 2-D subspace, and
     Steihaug CG.
2. `test_weights_bounds.f90`
   - vector versus matrix weights;
   - active lower parameter bound.
3. `test_robust.f90`
   - Huber IRLS on a contaminated straight line.
4. `test_multistart.f90`
   - BoxBOD quasi-random multi-start fit.
5. `test_large_operator.f90`
   - matrix-free Jacobian/Jacobian-transpose Steihaug-CG solve.
6. `test_stats_validation.f90`
   - hat values, Cook distance, log likelihood, traces, and invalid bounds.

All six tests pass.

## Numerical reference points

Misra1a converges near:

```text
par = 238.94212918, 5.5015643181e-4
RSS = 0.1245514
```

BoxBOD multi-start converges near:

```text
par = 213.80940889, 0.54723748542
RSS = 1168.0088766
```

A Huber fit to a line with one extreme outlier recovers parameters within
approximately `1e-5` of `(1,2)`, whereas ordinary least squares shifts the
intercept by roughly 0.81.

The matrix-free large-system regression recovers an exact two-parameter linear
solution with zero residual to floating-point precision.

## Build portability

The strict build is intentionally run with `-Werror=implicit-interface`. All
library callbacks are explicitly typed. Test/example programs use internal
procedures as callbacks; GNU/Linux may emit a linker note about an executable
stack because gfortran implements internal-procedure callbacks with trampolines.
That note is not an implicit-interface warning and does not affect the library
source.

FPM was not installed in the validation environment. `fpm.toml` is parsed with
Python's TOML parser, and the same source/test/example tree is compiled directly
by `scripts/test_gfortran.sh`.
