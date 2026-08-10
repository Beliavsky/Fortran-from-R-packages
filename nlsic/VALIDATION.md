# Validation

The translated source is validated with GNU Fortran using:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Permanent regression programs cover:

1. LDP feasibility/infeasibility and rank-deficient least-norm LSI.
2. Equality plus inequality constraints, custom norm matrices, null spaces,
   particular solutions, TLS, and numeric bound conversion.
3. The upstream constrained exponential nonlinear least-squares example with
   an analytical Jacobian.
4. Numerical-Jacobian nonlinear least squares with an equality constraint.
5. Rank-deficient regularized LS and infeasible nonlinear constraints.

The examples are also compiled and executed under the same flags.

Strict validation result for v0.1.0: **5/5 test programs passed** and both
examples compiled and ran successfully.  The constrained exponential example
returned `(1.0, 2.0)` with RSS approximately `5.8e-28`; the constrained
rank-deficient linear example returned `(0.25, 0.25, 0.5)`.

FPM itself was not installed in the validation container, so the same FPM
source tree was compiled directly with GNU Fortran.  `fpm.toml` was parsed
successfully with Python's TOML parser.
