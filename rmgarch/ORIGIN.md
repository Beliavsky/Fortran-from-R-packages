# Origin and provenance

- Original project: **rmgarch: Multivariate GARCH Models**
- Original package version: **1.4-2**
- Original package date: **2025-08-31**
- Original author and copyright holder: **Alexios Galanos**
- Original source archive supplied for this translation: `rmgarch-master(2).zip`
- Original archive SHA-256: `3db8063d8a9c4b69e20b2eefb4321b7bdf9d819f71d42ca7c9b9251ad9248197`
- Original package license declaration: **GPL-3**
- Translation license: **GPL-3.0-only**
- Translation status: **experimental and partially validated**

The original `DESCRIPTION` and `NAMESPACE` files are retained under
`reference/`. The supplied source archive did not expose a Git commit
identifier, so package version, package date, and archive hash identify it.

## Translation approach

The Fortran project independently reimplements selected numerical algorithms
in the R and native source. It does not embed R, Rcpp, Armadillo, `rugarch`, or
the original class system.

Plotting, S3/S4 methods, formula parsing, `xts`/`zoo` indexing, serialized R
datasets, parallel cluster management, and external package plumbing remain
outside the translation. Numerical omissions are listed explicitly in
`README.md`; no full-parity claim is made.
