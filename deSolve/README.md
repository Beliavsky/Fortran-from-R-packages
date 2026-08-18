# deSolve-fortran

A modern Fortran/FPM translation of the **computational core of R package
`deSolve` 1.42**.

The most important release constraint is explicit: **all Fortran source in
this package is free format**. The classic ODEPACK, VODE/ZVODE, DASPK, and
RADAU5 sources supplied by upstream deSolve were converted from fixed-column
Fortran into `.f90` free-form source. There are no `.f`, `.for`, or `.f77`
files and no fixed-form compiler option is required.

## Included solver families

The modern API exposes:

- `lsoda` -- ODEPACK automatic Adams/BDF switching
- `lsode` -- ODEPACK dense/banded ODE solver
- `lsodes` -- ODEPACK sparse solver
- `lsodar` -- LSODA with root detection and root-state recording
- `vode` and complex `zvode`
- `radau` -- Hairer/Wanner RADAU5
- `daspk` -- DAE solver
- `ode` -- simple method-name dispatcher
- generic Runge-Kutta integration with all 19 methods from upstream
  `rkMethod()` (including `ode23`/`ode45` aliases and the six experimental
  implicit methods)
- `euler` and `rk4` convenience entry points
- `dede_rk4` -- method-of-steps DDE driver using the upstream default Hermite
  history interpolation
- fixed-step discrete iteration via `iterate_map`

The original low-level solver entry points are also linked into the library,
so advanced callers can use ODEPACK work arrays/options directly if needed.

## Other translated computational utilities

- Brent root finding
- piecewise linear/constant forcing tables
- tabular events with replace/add/multiply semantics
- cubic Hermite lag-value and lag-derivative interpolation
- history buffer support
- 2-D and 3-D reaction/transport sparsity maps, returned as modern CSR
  patterns
- event-time cleaning and nearest-event utilities
- small self-contained BLAS-1 compatibility routines required by the
  classic solvers, avoiding an external BLAS dependency

## Build

```text
fpm build
fpm test
fpm run --example lotka_volterra
```

The manifest explicitly sets `source-form = "free"`. The converted classic
ODEPACK/VODE/DASPK/RADAU backends retain legacy implicit typing and implicit
external-procedure conventions, so the manifest also sets
`implicit-typing = true` and `implicit-external = true`. This is required by
current fpm, whose defaults disable both legacy features. The modern facade
modules themselves use explicit typing/interfaces.

A recent Fortran compiler supporting Fortran 2018 is recommended. The release
was validated with GNU Fortran 14.2.

## Minimal example

```fortran
use desolve, only : dp, ode_result, ode

type(ode_result) :: sol
real(dp) :: y0(1), times(3)

y0 = 1.0_dp
times = [0.0_dp, 0.5_dp, 1.0_dp]
sol = ode(rhs, y0, times, method='lsoda')
```

The right-hand side is a procedure with interface
`subroutine rhs(t, y, dydt)` using `real(dp)` arrays.

## Validation

The release tree was compiled with:

```text
-std=f2018 -O2 -fcheck=all
```

The modern API modules were additionally checked with:

```text
-Wall -Wextra -Wimplicit-interface -Werror
```

All **6/6 integration test programs pass**. Tests cover the ODEPACK/VODE
families, LSODAR roots, RADAU5, DASPK, ZVODE, all 19 Runge-Kutta tables,
DDE history, Brent roots, forcings/events, iteration, and 2-D/3-D sparsity
construction.

The `fpm` executable was not installed in the translation environment, so the
same source/test tree was compiled and executed directly with `gfortran`.
The `fpm.toml` file was parsed separately as TOML.

## License and provenance

The upstream package is GPL (>= 2); this translation is therefore distributed
as **GPL-2.0-or-later**. `COPYING`, upstream `DESCRIPTION`, `CITATION`, and
`NEWS.md` are retained. Converted numerical sources preserve the original
ODEPACK/VODE/DASPK/RADAU author and provenance comments.

See `TRANSLATION_NOTES.md` for the R/C/Fortran-to-Fortran coverage map and the
few API-level differences introduced by removing R's runtime/object system.
