# Licensing and provenance

## flexsurv

Source package: `flexsurv` 2.3.2.

Upstream `DESCRIPTION` states `License: GPL (>= 2)`.  The translated flexsurv
modules in this project are distributed under **GPL-2.0-or-later**.  The full
upstream source tree is retained at `upstream/flexsurv-2.3.2/`.

## numDeriv

The supplied `numDeriv-fortran` translation is derived from numDeriv
2016.8-1.1 (`GPL-2`).  The numerical derivative modules included in `src/`
retain their SPDX/license notices.  The supplied archive is preserved at
`vendor/numDeriv-fortran.zip`.

## deSolve

The supplied `deSolve-fortran-v0.1.1` translation is derived from deSolve 1.42
(`GPL (>= 2)`).  The Runge-Kutta modules compiled in `src/` retain their
SPDX/license notices.  The complete supplied archive is preserved under
`vendor/`.

## quadprog

The supplied `quadprog-fortran` translation is derived from quadprog 1.5-8
(`GPL (>= 2)`) and the Goldfarb-Idnani routines.  The included modules retain
their licensing headers and the supplied archive is preserved under `vendor/`.

## Other retained translations

The previously produced `mstate-fortran-v0.4.0` archive is retained under `vendor/` as an integration/provenance reference. Selected modules from the supplied `survival_f90` archive are compiled in v0.3.0 as described below. The relsurv rate-table exception is also described below.

## relsurv rate-table engine

Version 0.2.0 compiles `relsurv_kinds.f90` and `relsurv_ratetable.f90`, derived
from the translated computational rate-table layer of R package `relsurv` 2.3-3
(`GPL (>= 2)`).  These files retain GPL-2.0-or-later notices.  The complete
translation archive remains at `vendor/relsurv-fortran-v0.2.0.zip`.

## survival AFT modules

Version 0.3.0 compiles `survival_kinds.f90`, `survival_types.f90`,
`survival_linalg.f90`, and `survival_aft.f90` from the supplied survival Fortran
translation. These modules retain their **LGPL-2.0-or-later** notices and are used
to reproduce flexsurv's `survreg`-based starting-value path. The complete supplied
archive remains at `vendor/survival_f90.zip`.

## splines2-compatible natural spline

Version 0.3.0 adds `flexsurv_splines2ns.f90`, a native Fortran implementation of
the natural cubic B-spline construction used by `splines2::naturalSpline` for
flexsurv's optional `spline="splines2ns"` path. The compatibility module is
released under **GPL-3.0-or-later** and is attributed to the splines2 algorithmic
reference. Since GPL-3.0-or-later is a permitted later version of flexsurv's
GPL-2-or-later license, the combined v0.3.0 package is distributed under
**GPL-3.0-or-later**.
