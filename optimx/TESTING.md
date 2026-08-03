# Testing

Four tests cover:

1. Rosenbrock convergence with variable metric, conjugate gradient, and
   safeguarded Newton; bounded quadratic optimization.
2. Objective, gradient, Hessian, positive-definiteness, and KKT checks.
3. Multi-method optimization, multistart, polyalgorithm chaining, ranking, and
   result conversion.
4. Bounds adjustment, feasible step lengths, axial search, scaling checks, and
   solver-name discovery.

`test_gfortran.sh` builds with Fortran 2018 conformance, all common warnings as
errors, array bounds checking, backtraces, and floating-point traps. It then
rebuilds with `-O3 -Werror`.
