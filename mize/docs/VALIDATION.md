# Validation

The automated test program covers:

- Rosenbrock minimization with BFGS, SR1, and L-BFGS
- FR, PR+, HS+, DY, HZ+, and PR-FR nonlinear conjugate-gradient updates
- Exact-Hessian Newton and partial-Hessian optimization
- Hessian-vector truncated Newton
- Steepest descent, momentum, NAG, and delta-bar-delta
- Every accepted line-search name
- Stateful initialization and stepping
- Stored progress arrays
- Central finite-difference gradient diagnostics
- Polymorphic callback data and monitor cancellation

Build configurations used for release validation:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface -fcheck=all
-ffpe-trap=invalid,zero,overflow -fbacktrace
```

and

```text
-std=f2018 -O3 -Wall -Wextra -Wimplicit-interface
-Werror=implicit-interface
```

The optimized Rosenbrock example converges near `(0.999999752,
0.999999601)` with objective approximately `1.03e-12` in 33 iterations on GNU
Fortran 14.2.0. Exact iteration counts can vary with compiler and floating-point
settings.
