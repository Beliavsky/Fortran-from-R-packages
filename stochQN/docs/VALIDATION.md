# Validation

The automated test program exercises:

- the complete oLBFGS request cycle
- correction-pair creation and invalid memory-size handling
- oLBFGS convergence on a three-dimensional positive-definite quadratic
- SQN convergence using exact Hessian-vector products
- SQN convergence using large-batch gradient differences
- adaQN convergence using large-batch gradient differences and RMSProp scaling
- adaQN empirical-Fisher memory and correction-pair creation
- analytic logistic gradients against centered finite differences
- analytic logistic Hessian-vector products against gradient differences
- typed stochastic logistic fitting and prediction
- probability bounds and coefficient dimensions

For the independent quadratic problem

```text
A = diag(1, 2, 4)
b = (1, -2, 0.5)
```

the exact minimizer is `(1, -1, 0.125)`. The checked tests require oLBFGS to
reach a maximum absolute error below `2e-8`, and both SQN variants below
`1e-5`.

The logistic derivative tests use centered finite differences with step
`1e-6` and require maximum absolute discrepancies below `2e-9`.

Recommended compiler validation:

```text
gfortran -std=f2018 -Wall -Wextra -Wimplicit-interface \
  -Werror=implicit-interface -fcheck=all \
  -ffpe-trap=invalid,zero,overflow
```

## Original C-kernel comparison

The validation directory includes a small scalar BLAS compatibility header and
a driver that compiles the supplied C kernel directly. For five oLBFGS updates
on the identity quadratic, both the original C implementation and the Fortran
translation produce:

```text
x = (1.18098, -0.59049)
correction pairs = 4
```

The full double-precision values agree bit-for-bit in the tested environment.
Run `validation/c_reference/run_reference.sh` from the package root to reproduce
the C result. The corresponding Fortran check is part of `test/test_stochqn.f90`.
