# Reference generation

The test references were generated independently of the Fortran implementation.

- Fixed local-polynomial estimates and equivalent-kernel weights were computed in Python by constructing the weighted Vandermonde normal equations directly.
- Nadaraya-Watson values were computed from the package kernels by direct normalized weighted sums.
- The Buhlmann lag-window long-run variance reference was implemented independently from the recursion in the original Rcpp source.
- ARMA tests use seeded simulated series and check parameter recovery, residual identities, information criteria, and infinite-MA recursions rather than relying only on self-consistency.
- The hidden algorithm tables and polynomial kernels were recovered from the original package's `R/sysdata.rda` and recorded explicitly in the Fortran source.

The exact reference values and tolerances used by the executable suite are in `test/`.
