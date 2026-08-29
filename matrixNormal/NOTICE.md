# Notices and attribution

## matrixNormal

This is a modern Fortran translation of the computational code in `matrixNormal` 0.1.2 by Paul M. Hargarten.

Upstream metadata:

* Package: `matrixNormal`
* Version: 0.1.2
* Author/Maintainer: Paul M. Hargarten
* License: GPL-3
* Description reference: Pocuca, Gallaugher, Clark & McNicholas (2019), *Assessing and Visualizing Matrix Variate Normality*, arXiv:1910.02859.
* Distribution reference used by upstream documentation: Gupta & Nagar (1999), *Matrix Variate Distributions*.

The original attached source tree is retained unchanged under `upstream/matrixNormal-master`.

The matrix-square/symmetry/positive-definiteness helpers in upstream matrixNormal state that they were adapted from Frederick Novomestky's `matrixcalc` package. The `vec()` documentation likewise credits `matrixcalc` and cites Magnus & Neudecker (1999).

## mvtnorm Fortran dependency

The attached `mvtnorm-fortran` translation is retained under `vendor/mvtnorm-fortran`. Its own notices identify the upstream authors/contributors as Alan Genz, Frank Bretz, Tetsuhisa Miwa, Xuefei Mi, Friedrich Leisch, Fabian Scheipl, Bjoern Bornkamp, Martin Maechler, Torsten Hothorn, and contributors.

That dependency retains its own GPL-2.0-only license, source provenance, notices, and attribution. No mvtnorm source is relicensed by this port.

## Fortran port

Files in the top-level `src/` that derive from matrixNormal are marked `SPDX-License-Identifier: GPL-3.0-only`.
