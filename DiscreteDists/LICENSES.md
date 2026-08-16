# Licensing

## DiscreteDists-derived Fortran source

The supplied `DiscreteDists` package declares `MIT + file LICENSE`.  The new
Fortran translation in `src/`, tests, examples and accompanying translation
documentation are distributed under the same MIT terms.  See `LICENSE`.

## Vendored COMPoissonReg-fortran dependency

`vendor/compoissonreg-fortran` is a separately identified dependency translated
from `COMPoissonReg`, whose upstream license is `GPL-2 | GPL-3`.  Its own
license files and provenance documentation remain inside that directory.
Programs linked with that dependency must comply with the applicable GPL terms.

## Upstream source snapshot

The complete supplied `DiscreteDists` source is retained under `upstream/`
with its original metadata and license files unchanged.

The supplied `gamlss-fortran`, `gamlss.dist-fortran`, and `nleqslv-fortran`
archives are not copied into this release and are not runtime dependencies.
