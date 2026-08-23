# Dependency map

The upstream R package imports `bigstatsr`, `Rfast`, `Rfast2`, `Rnanoflann`, and `rangen` among GUI/geospatial/parallel packages.
The supplied Fortran translations were inspected during this port.

- **Rfast-fortran v0.3.0** already contains `rfast_directional` (`vm_mle`, `vmf_mle`, `rvonmises`, Watson/Kuiper, circular-linear correlation), distributions, distances, regression, matrix utilities and MLE infrastructure. These are natural backends for future exact parity of Directional wrappers.
- **Rfast2-fortran v0.1.0** contains random sampling, column MLE, multivariate and regression/statistics modules corresponding to many `Rfast2::` calls.
- **Rnanoflann-fortran v0.1.0** provides nearest-neighbor search primitives suitable for `dirknn`, `cosnn`, tuning and KNN regression parity.
- **bigstatsr-fortran v0.1.0** provides FBM/matrix/statistics/SVD/regression facilities corresponding to `read.fbm` and large-data pathways.
- **rangen-fortran v0.1.0** provides PCG32 and distribution/sampling routines useful for reproducible bootstrap, permutation and simulation routines.

v0.1.0 intentionally does not vendor or link these packages, avoiding duplicate kind/module/runtime assumptions and keeping this package independently buildable. The numerical formulas directly translated here are derived from Directional itself.
