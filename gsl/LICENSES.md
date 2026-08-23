# Licensing and attribution

## gsl-fortran translated interface

SPDX-License-Identifier: GPL-3.0-only

This project translates the computational interface of the R package `gsl` 2.1-9. The upstream package `DESCRIPTION` declares `License: GPL-3`. The GNU GPL version 3 text is included as `LICENSE`.

## Upstream R package

- Package: `gsl`
- Version: 2.1-9
- Author/maintainer: Robin K. S. Hankin
- Contributors named upstream: Andrew Clausen and Duncan Murdoch
- License: GPL-3
- System requirement: GNU Scientific Library >= 2.5

The supplied source snapshot is preserved unchanged in `upstream/gsl-master.zip`.

The C special-function shims in `src/c_shim/` are derived directly from the upstream package's thin C wrappers and therefore remain under the same GPL-3 terms. `fortran_rng.c` and `fortran_qrng.c` are new interoperability shims implementing the same GSL object operations without R's `SEXP` layer.

## GNU Scientific Library

GNU GSL is not vendored in this archive. It is an external linked dependency. Users must install a compatible GSL development package. Its licensing remains governed by the GNU GSL distribution.
