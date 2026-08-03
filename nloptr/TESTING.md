# Testing

Four test programs cover:

1. projected BFGS, Nelder-Mead, and derivative-free wrapper behavior on the
   Rosenbrock function;
2. nonlinear inequalities, nonlinear equality handling, and bound projection;
3. central finite-difference gradients/Jacobians and derivative checking;
4. deterministic global multistart, options helpers, and result validation.

Validation commands:

```text
./scripts/test_gfortran.sh
./scripts/test_gfortran_optimized.sh
```

The strict configuration uses Fortran 2018, bounds checking, floating-point
traps, and warnings as errors. The optimized configuration uses `-O3 -Werror`.
