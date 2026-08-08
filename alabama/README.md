# alabama-fortran

Modern Fortran 2018 computational port of R package **alabama 2025.1.0**.

The port implements both constrained optimization families in the package:

- `auglag`: augmented-Lagrangian method supporting equality constraints,
  inequality constraints, or both, and accepting infeasible starting points.
- `constr_optim_nl`: the older adaptive-barrier / augmented-penalty method.
- Public helper routines corresponding to the R implementation:
  `auglag1`, `auglag2`, `auglag3`, `adpbar`, `augpen`, and `alabama_legacy`.

Inequalities follow the R convention `hin(x) >= 0`; equalities use `heq(x) = 0`.
The legacy adaptive-barrier routines require `hin(x) > 0` at the initial point.

## Numerical dependencies

The supplied `numDeriv-fortran` port is bundled as an FPM path dependency and is
used when objective gradients or constraint Jacobians are not supplied. It is
also used for the augmented-Lagrangian Hessian/KKT check.

The R package delegates inner minimization to `stats::optim` or `nlminb`. This
Fortran port uses the bundled `roptim` translation to provide BFGS,
Nelder-Mead, nonlinear CG, L-BFGS-B and SANN. The string `"nlminb"` is accepted
for source compatibility and mapped to BFGS because R's `nlminb` implementation
is not part of the alabama source package.

## Build

```text
fpm build
fpm test
fpm run --example constrained_example
```

The package name in `fpm.toml` is `alabama`. Dependency names exactly match the
bundled dependency manifests, avoiding FPM dependency-name mismatches.

Version 0.1.1 also routes all user callbacks through module-level procedures with
explicit abstract interfaces. This avoids GNU Fortran/FPM `implicit-interface`
diagnostics for optional callbacks captured by internal optimizer procedures.

## Minimal example

```fortran
use alabama
real(dp) :: x(2)
type(alabama_result_t) :: fit

x = [-0.2_dp, 2.6_dp]
call auglag3(x, objective, inequalities, equality, fit)
```

See `example/constrained_example.f90` for complete callbacks.

## Translation scope

R lists, S3-like result printing, `...` argument forwarding, and calls to the R
runtime are intentionally omitted. Numerical controls and result fields are
represented by typed Fortran derived types. Plotting is not part of alabama.
