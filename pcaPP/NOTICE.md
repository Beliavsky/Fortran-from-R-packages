# Notices and provenance

This directory is a clean-room-style Fortran re-expression of the
computational behavior of **pcaPP 2.0-5**, translated from the source archive
supplied with this task.  The upstream package identifies Peter Filzmoser,
Heinrich Fritz, Klaudius Kalcher, and Valentin Todorov as authors and is
licensed under GPL (>= 3).  The Fortran translation is distributed under
GPL-3.0-or-later; see `LICENSE`.

The files under `upstream/` are retained only to document the exact R package
revision and the R functions used as translation references.  The original
C/C++ implementation is not copied or vendored into this package.

The upstream native implementation contains substantial pcaPP code developed
by Heinrich Fritz, including robust PCA and L1-median routines.  Those
algorithms were independently re-expressed here in modern free-form Fortran.
The Fortran implementation deliberately uses typed numeric APIs instead of R
S3 objects and `.Call`/`.C` interfaces.

The upstream `src/cov.kendall.cpp` contains David Simcha's 2010 Kendall tau
implementation under the Boost Software License 1.0 and cites W. R. Knight,
"A Computer Method for Calculating Kendall's Tau with Ungrouped Data",
JASA 61 (1966), 436-439.  This translation does not copy that source: it uses
an independent O(n^2) tau-b implementation.  Attribution and the Boost 1.0
license are nevertheless preserved here because that code is part of the
upstream provenance; see `LICENSE.BOOST`.

`qn` follows pcaPP's normalization and finite-sample correction conventions,
but computes the required pairwise-distance order statistic directly rather
than copying the upstream selection code.

No BLAS, LAPACK, R runtime, `mvtnorm`, `robustbase`, or other translated R
package source is embedded in this package.  The current numerical surface is
self-contained.
