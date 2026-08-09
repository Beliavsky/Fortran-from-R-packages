# quadprogXT-fortran

Modern Fortran/FPM translation of the computational code in the R package
`quadprogXT` 0.0.6.

`quadprogXT` extends ordinary convex quadratic programming with linear
constraints and linear penalties involving absolute values of the decision
variables and their changes from a reference point.  The translation performs
the same auxiliary-variable transformation and solves the resulting QP using
the supplied `quadprog-fortran` Goldfarb-Idnani implementation.

## Features

- `solve_qp_xt`: build and solve the extended QP in one call.
- `build_qp_xt`: expose the transformed ordinary QP.
- `normalize_constraints`: normalize each constraint column and its bound by
  the column 2-norm.
- `convert_to_compact`: convert a dense constraint matrix to the compact
  representation consumed by `solve_qp_compact`.
- Absolute-value constraints/penalties on `b`.
- Absolute-change constraints/penalties on `b-b0`.
- Dense or compact constraint solving.
- Constraint normalization on/off.
- Equality constraints through `meq`.
- Factorized `Dmat` support matching `quadprog`.

## Build

```text
fpm build
fpm test
```

The supplied `quadprog-fortran` translation is vendored as an FPM path
dependency under `vendor/quadprog-fortran`.

## Example: L1 constraint

```fortran
use quadprog_kinds, only: dp
use quadprog, only: qp_result
use quadprogxt, only: solve_qp_xt

real(dp) :: dmat(2,2), dvec(2), aabs(4,1), babs(1)
type(qp_result) :: fit

dmat = 0.0_dp
dmat(1,1) = 1.0_dp
dmat(2,2) = 1.0_dp
dvec = [2.0_dp, 1.0_dp]
aabs(:,1) = -1.0_dp
babs = -1.0_dp

fit = solve_qp_xt(dmat, dvec, amat_posneg=aabs, bvec_posneg=babs)
```

The first two elements of `fit%solution` are the original decision vector.
Additional elements are the auxiliary absolute-value variables introduced by
`quadprogXT`.

See `API.md` and `TRANSLATION_COVERAGE.md` for details.
