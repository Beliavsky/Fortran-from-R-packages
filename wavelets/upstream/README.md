# Retained upstream computational sources

This directory contains source snapshots from R package `wavelets` 0.3-0.2 used
as the provenance basis for the modern Fortran translation.

`DESCRIPTION` and `NAMESPACE` record the upstream metadata/export policy.
`R/` contains only the counted computational R entry points used for translation.
`src/` contains the four upstream C transform kernels translated into pure
Fortran. Plotting, printing, summary, and other presentation/interface sources
are intentionally not copied into this translation package.

These retained files remain subject to the upstream package's GPL (>= 2)
license declaration. See the repository-level `COPYING`, `LICENSE`, and `NOTICE`.
