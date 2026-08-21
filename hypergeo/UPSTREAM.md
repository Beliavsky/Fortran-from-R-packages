# Upstream provenance

This project translates the computational code of:

- R package: `hypergeo`
- Version: 1.2-14
- Title: The Gauss Hypergeometric Function
- Author: Robin K. S. Hankin
- Upstream license: GPL-2
- Upstream package description: Gaussian hypergeometric function for complex numbers

The original `DESCRIPTION`, `NAMESPACE`, R computational sources, tests, manual pages,
and vignette sources are retained under `upstream/` for provenance and comparison.

The translation also vendors previously translated dependencies:

- `elliptic-fortran` 0.1.0, corresponding to the R `elliptic` dependency
- `contfrac-fortran` 0.1.0, corresponding to the R `contfrac` dependency
- `desolve` 0.1.1, corresponding to the R `deSolve` dependency

Each vendored dependency retains its own license/provenance files.
