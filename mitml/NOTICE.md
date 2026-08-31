# NOTICE

## Upstream package

This work is derived from the computational routines of R package **mitml
0.4-5**, dated 2023-03-08.

Upstream metadata names:

- Simon Grund -- author and creator/maintainer
- Alexander Robitzsch -- author
- Oliver Luedtke -- author

The upstream `DESCRIPTION` declares `License: GPL (>= 2)`. The original
`DESCRIPTION`, `NAMESPACE`, `README.md`, and the R files used as direct
computational references are retained byte-for-byte under `upstream/`.

No copyright year or copyright holder beyond the attribution present in the
upstream metadata has been invented for this translation.

## Translation

The maintained Fortran source is distributed under GPL-2.0-or-later, consistent
with the upstream license. The full GNU GPL version 2 text is included in
`LICENSE`; the upstream grant permits version 2 or any later GPL version.

The translation replaces R-specific model/list/formula infrastructure with a
numeric array API. It does not contain source copied from `pan`, `jomo`,
`rfortran-core`, `rfortran-linalg`, BLAS, or LAPACK.

`pan` and `jomo` are documented companion translations rather than simultaneous
FPM dependencies. This avoids creating a single statically linked dependency
graph combining the current GPL-3-only `pan` translation and GPL-2-only `jomo`
translation.
