# Third-party and upstream notices

This project is a modern-Fortran translation of computational ideas and
algorithms from the R package `LaplacesDemon` 16.1.8.

The upstream package is distributed under the MIT license. Its license and
package metadata are retained with the translation, and selected original R
files are kept under `upstream/reference-R*` for algorithm/provenance
comparison.

The upstream package notes that its trust-region (`TR`) method in
`LaplaceApproximation` is derived from `trust::trust` by Charles J. Geyer. The
Fortran `trust_region_maximize` routine is an independent numerical
implementation of a standard trust-region Newton strategy and does not copy
`trust::trust` source code. The retained upstream R wrapper is present only for
provenance/comparison under the upstream MIT licensing terms.

No R runtime, BLAS-specific extension, C/C++ bridge, or other third-party
binary library is required by this release.
