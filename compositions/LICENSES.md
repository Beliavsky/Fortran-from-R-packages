# Licenses and provenance

## compositions

The translated package is derived from `compositions` 2.0-9, licensed GPL
(version 2 or later).  The complete supplied upstream tree is retained under
`upstream/`, including its `COPYING` file.

New translation modules in `src/compositions_*.f90` are marked
`SPDX-License-Identifier: GPL-2.0-or-later`.

## robustbase translation

The supplied `robustbase-modern-fortran-0.3.0` translation is retained under
`vendor/robustbase-modern-fortran-0.3.0`.  The subset compiled into `src/`
retains its original SPDX/license headers.

## bayesm translation

The supplied `bayesm-fortran-v0.1.0` translation is retained under
`vendor/bayesm-fortran-v0.1.0`.  The subset compiled into `src/` retains its
original SPDX/license headers.

## tensorA translation

The supplied `tensorA-fortran-v0.1.0` translation is retained under
`vendor/tensorA-fortran-v0.1.0`. Its Fortran modules are also compiled into the
main library for named high-rank tensor operations; their original GPL headers
are retained.
