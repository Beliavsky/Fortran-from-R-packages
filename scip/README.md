# scip-fortran

Modern Fortran/FPM translation of the computational interface in R package
**scip 1.10.0-3**, using the package's exact vendored **SCIP 10.0.2** and
**SoPlex 8.0.2** sources as the optimization backend.

## What is translated

The Fortran API covers the full exported computational surface of the R
package:

- one-shot LP/MIP solving from dense or CSC matrices;
- continuous, binary, and integer variables;
- persistent incremental model construction;
- ranged linear constraints;
- quadratic constraints;
- SOS1 and SOS2 constraints;
- indicator constraints;
- objective minimization/maximization;
- common SCIP limits, tolerances, presolve, LP, branching, heuristic, and
  thread controls;
- arbitrary additional native SCIP parameters;
- solver status, objective, primal solution, solution-pool access, node count,
  LP iterations, solve time, gap, and solution count.

R S3 printing, `Matrix`/`slam` class adapters, and `.Call` marshalling are
omitted.  There is no plotting code in the package's solver layer.

## Why SCIP itself is still C/C++

The uploaded R package vendors the complete SCIP Optimization Suite.  SCIP is a
large independent solver framework, not package-specific wrapper code.  This
project translates the R-facing computational API into Fortran and uses
`iso_c_binding` plus a small plain-C shim to call the **same bundled solver**.
It does not replace SCIP with a smaller LP/MIP implementation, and it does not
claim that SCIP's hundreds of thousands of lines have been mechanically
rewritten into Fortran.

## Build

The final archive is source-only.  Build the bundled solver once, then use FPM.

### Linux / macOS

```sh
./scripts/build_vendor.sh
./scripts/fpm_build.sh
./scripts/fpm_test.sh
```

`JOBS=8 ./scripts/build_vendor.sh` can be used for a parallel backend build.
For a quicker non-optimized backend build while developing, set
`SCIP_FORTRAN_FAST_BUILD=1`.

### Windows with MinGW/gfortran

From PowerShell, with CMake, GNU `gcc`/`g++`/`ar`, and FPM on `PATH`:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\build_vendor.ps1
powershell -ExecutionPolicy Bypass -File scripts\fpm_build.ps1
```

The PowerShell backend builder deliberately uses the **MinGW Makefiles** CMake
generator so the C/C++ static library is ABI-compatible with gfortran.

### Plain `fpm build`

After `vendor/lib/libscipfortran_backend.a` has been built, `fpm.toml` links it
as `scipfortran_backend`.  If invoking FPM directly, ensure `vendor/lib` is on
GCC's `LIBRARY_PATH`.  The supplied `fpm_build` scripts do this automatically.

## One-shot example

```fortran
use scip
implicit none
real(dp) :: a(2,2), obj(2), b(2)
character(len=2) :: sense(2)
type(scip_control) :: ctrl
type(scip_result) :: res

a = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2,2])
obj = [-5.0_dp, -4.0_dp]
b = [6.0_dp, 8.0_dp]
sense = ['<=', '<=']
ctrl%verbose = .false.

res = scip_solve(obj, a, b, sense, control=ctrl)
print *, trim(res%status)
print *, res%objval
print *, res%x
```

This gives objective `-22` at `(10/3, 4/3)`.  As in the R package,
`scip_solve` minimizes; negate the objective for maximization.

## Model-building example

```fortran
use scip
implicit none
type(scip_model_t) :: m
type(scip_solution) :: sol
integer :: first, cidx

m = scip_model('choose_one')
call m%set_objective_sense('maximize')
first = m%add_vars([3.0_dp,2.0_dp,1.0_dp], &
                   [0.0_dp,0.0_dp,0.0_dp], &
                   [1.0_dp,1.0_dp,1.0_dp], ['B','B','B'])
cidx = m%add_linear_cons([1,2,3], [1.0_dp,1.0_dp,1.0_dp], &
                          lhs=1.0_dp, rhs=1.0_dp)
call m%optimize()
sol = m%get_solution()
print *, sol%objval, sol%x
call m%free()
```

`scip_model_t` is an opaque owning handle.  Do not copy it with intrinsic
assignment after construction, and explicitly call `%free()` when finished.

## CSC matrices

`type(scip_csc_matrix)` stores 1-based Fortran row indices and 1-based column
pointers.  `make_csc_matrix(dense)` converts a dense matrix.  The generic
`scip_solve` accepts either a dense rank-2 array or this CSC type.

## Tests

Five regression programs cover dense/CSC LP, binary MIP, incremental models,
solution-pool queries, quadratic constraints, SOS1/SOS2, indicators,
infeasible/unbounded statuses, and controls.  See `VALIDATION.md`.

## License

Apache-2.0.  See `LICENSE` and `LICENSES.md`.  The complete uploaded R package,
including the original SCIP and SoPlex source/license files, is retained under
`original/`.
