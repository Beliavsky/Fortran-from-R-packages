# Upstream provenance

This project translates the computational code of R package `rootSolve` 1.8.2.4.

Upstream package metadata identifies:

- Karline Soetaert — author and maintainer
- Alan C. Hindmarsh — LSODES/sparse contributions
- S. C. Eisenstat — sparse solver contribution
- Cleve Moler and Jack Dongarra — LINPACK contribution
- Youcef Saad — SparseKit contribution

The upstream package declares `License: GPL (>= 2)`.

The release retains the upstream `DESCRIPTION`, `NAMESPACE`, R computational sources, citation, and R interface C sources under `upstream/rootSolve-master/`.  The original fixed-form Fortran files are not duplicated there because this translation has an explicit all-free-format Fortran requirement.  The mathematical algorithms and author provenance of those files are documented in `TRANSLATION_NOTES.md` and in the translated/vendored free-format sources.

The ODEPACK backend used by `runsteady` comes from the previously translated free-format `deSolve` package, itself GPL-2.0-or-later and retaining the original solver attribution comments.
