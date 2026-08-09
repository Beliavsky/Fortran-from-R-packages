# lpSolve-fortran

Modern Fortran/FPM translation of the exported computational core of the R
package `lpSolve` 5.6.23.9000 / lp_solve 5.5.

## Implemented

- General continuous LP minimization and maximization.
- Mixed-integer and binary LPs by branch-and-bound.
- Multiple all-binary solutions using exact no-good constraints.
- `<=`, `>=`, and `=` constraints.
- Dense and sparse-triplet constraint input.
- Nonnegative decision variables, matching the R front end.
- Assignment-problem helper corresponding to `lp.assign`.
- Transportation-problem helper corresponding to `lp.transport`.
- `make.q8` 8-queens sparse constraint generator.
- LP row duals and reduced costs.
- Row equilibration, iteration/node limits, and timeout controls.
- lp_solve-compatible primary status codes for optimal, suboptimal,
  infeasible, unbounded, numerical failure, and timeout.

The default numerical core is a standalone two-phase primal simplex tableau.
MILPs use LP-relaxation branch-and-bound.  No C, BLAS, or LAPACK dependency is
required.

## Build

```text
fpm build
fpm test
fpm run --example basic_lp
fpm run --example assignment
```

## Example

```fortran
use lpsolve

type(lp_result) :: result
real(dp) :: c(3), a(2,3), b(2)
integer :: sense(2)

c = [1.0_dp, 9.0_dp, 1.0_dp]
a = reshape([1.0_dp, 3.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 2.0_dp], [2,3])
b = [9.0_dp, 15.0_dp]
sense = LP_LE

call solve_lp(LP_MAX, c, a, sense, b, result)
```

For the upstream documented example this returns objective `40.5` and solution
`[0, 4.5, 0]`.  Marking all variables integer returns objective `37` and
solution `[1, 4, 0]`.

See `TRANSLATION_COVERAGE.md` for the advanced lp_solve subsystems that are not
part of this first standalone Fortran release.
