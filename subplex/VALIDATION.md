# Validation

The release was validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression coverage includes:

1. 2-D and 6-D Rosenbrock tests from upstream.
2. 1-D, 2-D, and 3-D ripple tests from upstream.
3. Shifted Rosenbrock solution and numerical Hessian from upstream.
4. Maximum-evaluation, loose-tolerance, scalar-scale, and vector-scale cases.
5. A direct numerical-Hessian test on a known quadratic.

Both examples are also compiled and executed by the strict validation script.

For an additional source-fidelity check, the original fixed-form Rowan routines
were compiled directly with GNU Fortran and BLAS. The 2-D Rosenbrock optimum,
objective value, evaluation count, and convergence code matched the modern
translation exactly in this environment. The 3-D ripple case reached the same
published local minimizer to the upstream test tolerance; its exact evaluation
trajectory differed slightly because the modern translation replaces BLAS
copy/scale/AXPY operations and COMMON storage with intrinsic array arithmetic.

FPM itself was not installed in the validation container. `fpm.toml` was parsed
with a TOML parser, and the exact FPM source/test/example tree was compiled by
the strict scripts.
