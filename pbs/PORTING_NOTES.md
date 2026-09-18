# Porting notes

## Numerical design

`pbs` uses a Cox-de Boor recurrence for normalized B-spline values. Derivatives needed by the nonperiodic extrapolation path are generated recursively from lower-degree basis derivatives, matching the Taylor expansion used by upstream `pbs()` outside explicit ordinary-spline boundary knots.

For periodic splines, the boundary/internal knot vector is extended by `degree` knots on each side using the same wrapped spacing as upstream. The first and last `degree` basis columns are then added pairwise, yielding a basis whose values at the two period boundaries coincide.

When `df` is supplied and explicit knots are absent, knot probabilities follow the upstream formulas and quantiles use R's default type-7 interpolation rule.

## Interface choices

Fortran does not have R's S3 attributes or lazy formula objects. `type(pbs_basis)` therefore stores the matrix and reusable spline specification explicitly. `predict_pbs` is the computational counterpart of `predict.pbs`.

`makepredictcall.pbs` has no numerical analogue and is intentionally not translated.
