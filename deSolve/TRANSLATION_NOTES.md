# Translation notes: deSolve 1.42 -> deSolve-fortran 0.1.0

## Scope

This project translates the computational code of the attached `deSolve`
1.42 source tree to a standalone FPM library. R object management, S3 methods,
plotting, dynamic-library registration, and R-specific `.Call`/`.C` glue are
not part of the numerical library.

## Free-format requirement

Every Fortran source shipped in this release is free format (`.f90`). In
particular, the following upstream fixed-format files were converted to valid
free-form syntax while retaining their numerical statements and comments:

- `daux.f`
- `ddaspk.f`
- `dintdy2.f`
- `dlinpk.f`
- `dlsoder.f`
- `dsparsk.f`
- `dvode.f`
- `errmsg.f`
- `opkda1.f`
- `opkdmain.f`
- `radau5.f`
- `radau5a.f`
- `zvode.f`

The conversion handles fixed-column statement labels, continuation columns,
comment cards, and fixed-form blank-insensitive dotted operators. The result
compiles as free form with `-std=f2018`; `-ffixed-form` is neither used nor
needed.

## Classic solver backends

The numerical arithmetic of ODEPACK, VODE/ZVODE, DASPK, and RADAU5 is kept
close to the deSolve-vendored source because these are mature, extensively
validated algorithms. Modern Fortran modules provide procedure interfaces,
derived result types, automatic workspace allocation, and status/statistics
extraction around those backends.

R-specific printing/error calls introduced by deSolve were replaced with
standalone Fortran `error_unit`/`error stop` helpers. The level-1 BLAS symbols
used by the classic solvers are supplied by small reference Fortran routines,
so the FPM package has no external BLAS dependency.

## Modern high-level solver API

- `desolve_odepack`: `lsoda`, `lsode`, `vode`
- `desolve_roots_sparse`: `lsodes`, `lsodar`
- `desolve_stiff`: `radau`, `daspk`, `zvode`
- `desolve_rk`: generic Runge-Kutta solver and all tables from `rkMethod.R`
- `desolve_dde`: Hermite-history method-of-steps DDE solver
- `desolve`: public facade and `ode` method dispatcher

The high-level ODEPACK/VODE wrappers currently request internally generated
Jacobians. R's ability to pass arbitrary R/DLL Jacobian, mass-matrix,
preconditioner, or sparse-structure callbacks is not reproduced by those
convenience wrappers. The converted low-level solver entry points are still
present for callers that need the complete original option surface.

`DLSODER` and `DLSODESR`, which are internal/non-exported variants in deSolve,
are converted and built as part of the backend even though no separate modern
wrapper is provided.

## Runge-Kutta C code

The computational behavior of `rk_fixed.c`, `rk_auto.c`, and
`rk_implicit.c` was translated to native modern Fortran rather than retaining
C. `rk_method_by_name()` contains all 19 upstream method tables:

`euler`, `rk2`, `rk4`, `rk23`, `rk23bs`, `rk34f`, `rk45f`, `rk45ck`,
`rk45e`, `rk45dp6`, `rk45dp7`, `rk78dp`, `rk78f`, `irk3r`, `irk5r`,
`irk4hh`, `irk6kb`, `irk4l`, and `irk6l`.

`ode23` aliases `rk23bs`; `ode45` aliases `rk45dp7`, as upstream does.
Adaptive methods use the embedded pair and the same basic safety/min/max step
scaling constants. Because the Fortran API steps exactly to requested output
times, the R/C dense-output buffering machinery is unnecessary for ordinary
`rk_integrate` calls.

The implicit RK implementation solves the coupled stage equations by Newton
iteration with a finite-difference Jacobian and dense pivoted elimination,
matching the mathematical method while replacing the C/LINPACK glue.

## Forcings and events

The computational portions of `forcings.c` are represented by:

- `forcing_table`: linear or piecewise-constant interpolation
- `event_table`: replace, add, and multiply event semantics
- `clean_event_times` and `nearest_event`

R's event data-frame/list parsing and compiled-DLL function pointers are R
runtime infrastructure and are omitted. Event application is explicit in the
Fortran API rather than automatically mutating an R solver call.

## Delays

The package-default delay interpolation (`interpol=1`) is the cubic Hermite
algorithm from `lags.c`; `history_buffer%lag_value` and `%lag_deriv` implement
those formulas. The converted backends also retain `dintdy2` and RADAU dense
interpolation support used by deSolve's higher-order history modes.

The modern `dede_rk4` driver is a method-of-steps RK4 implementation. It is
not a wrapper that injects history updates into every internal LSODA/VODE
accepted step, because that mechanism in upstream deSolve depends on global C
state and R callbacks. For positive delays, callers should choose a step size
small relative to the shortest delay; extrapolation from the latest Hermite
history point is used when a stage requests a time newer than the last
accepted point.

## Sparse PDE maps

`twoDmap.c` encoded neighbor structure directly into ODEPACK's `IWORK` layout.
The translation exposes the same reaction/nearest-neighbor/cyclic-boundary
connectivity as `sparsity_2d` and `sparsity_3d`, returned as a conventional
CSR `sparsity_pattern`. This is more useful in standalone Fortran while the
low-level DLSODES interface remains available for users needing its exact
workspace convention.

## Discrete iteration and Brent root finding

`call_iteration.c` becomes `iterate_map`; the R model callback is replaced by
a typed Fortran procedure. `brent.c` becomes `brent_root` using the same
Brent/Forsythe-Malcolm-Moler iteration.

## Omitted R-only layers

The following are intentionally not translated as numerical Fortran code:

- S3 `print`, `plot`, `image`, `hist`, `summary`, `subset`, and diagnostics
  presentation methods
- R list/formula/name/dimension reshaping and argument validation machinery
- DLL registration/external-pointer code (`R_init_deSolve.c`, `DLLutil.c`)
- R-to-C callback marshaling (`call_*.c`)
- plotting (`matplot.R` and plotting helpers)
- bundled demonstration model DLLs, which are examples of the old R native
  extension interface rather than solver algorithms

The example and tests in this project are native Fortran replacements.


## v0.1.1 FPM compatibility fix

Current fpm defaults to disabling implicit typing and implicit external
interfaces. The converted classic solver backends intentionally preserve these
legacy Fortran conventions even though every source file is free format. The
manifest now explicitly sets:

```toml
[fortran]
source-form = "free"
implicit-typing = true
implicit-external = true
```

DASPK's legacy workspace-index `PARAMETER` constants were additionally changed
to explicit `integer, parameter ::` declarations. No fixed-form source was
reintroduced.
