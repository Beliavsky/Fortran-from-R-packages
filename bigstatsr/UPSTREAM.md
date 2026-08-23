# Upstream sources

This port was prepared from the supplied CRAN source tree:

- Package: `bigstatsr`
- Version: 1.6.2
- Date: 2025-06-30
- License: GPL-3
- Author/maintainer: Florian Prive, with Michael Blum and Hugues Aschard

The exact supplied archive is stored as `upstream/bigstatsr-master.zip`.

`bigstatsr` uses RSpectra for partial eigendecomposition/SVD. The Fortran port
vendors the previously translated RSpectra computational layer under MPL-2.0;
the corresponding source archive is retained as
`upstream/RSpectra-fortran-v0.1.0.zip`.
