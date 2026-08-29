# rfortran-compat

This transitional MIT-licensed package provides the compatibility runtime used by translations generated with R-to-Fortran. It allows older translations to share one runtime while its generally useful procedures are progressively migrated into the focused `rfortran-*` modules.

New translation code should prefer the focused shared modules. This package preserves the generated API for existing translations.
