# Validation

Validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. finite-difference gradient, Hessian, Jacobian, and Richardson gradient;
2. Rosenbrock convergence for BFGS, L-BFGS-B, Newton-Raphson, modified Newton,
   dogleg, double-dogleg, and Levenberg-Marquardt;
3. Gauss-Newton and least-squares modes of LM/dogleg/double-dogleg;
4. active box constraints for L-BFGS-B;
5. rejection of a stationary saddle when positive-definiteness checking is on;
6. explicit `gn_hessian` curvature paths for LM and double-dogleg.

The GNU linker may print an executable-stack note for tests/examples that pass
contained procedures as callbacks. This is produced by gfortran's trampoline
implementation for internal procedures; it is not an implicit-interface warning
and does not occur when callbacks are ordinary module procedures.
