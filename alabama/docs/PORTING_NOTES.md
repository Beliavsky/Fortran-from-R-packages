# Porting notes

## Derivatives

When `gr`, `hin_jac`, or `heq_jac` are omitted, the port calls the bundled
`numDeriv-fortran` simple finite-difference routines, matching the R package's
explicit `method="simple"` fallbacks. `alabama_inner_control_t%ndeps` controls
the finite-difference step.

## Inner optimizer

R `alabama` calls `stats::optim` (and optionally `nlminb`). The source package
does not contain those algorithms. This port uses the separately licensed
`roptim` Fortran package for the `optim`-style methods. `nlminb` is accepted as
a method name but maps to BFGS; this is documented rather than presented as a
direct `nlminb` translation.

## Maximization

R uses negative `fnscale`. The Fortran control exposes `maximize=.true.` and
applies the same sign convention to both the inner optimizer and penalty terms.

## Scaling

`i_scale` and `e_scale` are allocatable vectors. An unallocated component means
unit scaling; a length-one vector is broadcast; otherwise the vector length must
match the corresponding number of constraints.

## KKT2

The R package computes eigenvalues of the augmented-Lagrangian Hessian. The
Fortran port symmetrizes the numerically differentiated Hessian and tests
positive definiteness (or negative definiteness for maximization) by Cholesky
factorization, which is equivalent for a real symmetric nonsingular Hessian.

## Fortran-specific robustness correction

Two R expressions use an `if (is.infinite(sig0))` guard after forming
`sig0 = sig/Kprev`. In Fortran, an eager `merge(sig/Kprev,...)` can still
evaluate the zero denominator. The port uses an explicit branch when `Kprev=0`
so builds with floating-point traps are safe.


## Callback portability (0.1.1)

Some GNU Fortran/FPM configurations diagnose a direct call to an optional procedure
dummy from a host-associated internal procedure as an implicit-interface call. The
0.1.1 port therefore invokes objectives, gradients, constraints, and constraint
Jacobians through module-level trampoline routines whose dummy procedures are
declared with the package abstract interfaces. This changes no numerical behavior
or public callback signatures.
