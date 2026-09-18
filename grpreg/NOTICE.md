# NOTICE and provenance

This directory is a modern Fortran translation of computational functionality from the R package **grpreg 3.6.0**, "Regularization Paths for Regression Models with Grouped Covariates".

Upstream authors listed by the package are Patrick Breheny, Ryan Miller, Yaohui Zeng, and Ryan Kurth. The upstream package is distributed under GPL-3. The complete GPL version 3 text is retained in `LICENSE`.

The upstream package metadata, R computational sources, native C sources, CITATION, NEWS, and README used during translation are retained under `upstream/` so that numerical provenance and attribution can be inspected without mixing the original implementation into the FPM build.

The Fortran translation reimplements the group-descent and local-coordinate-descent algorithms from the upstream C sources, including group lasso, group MCP, group SCAD, group exponential lasso, composite MCP, group bridge, Gaussian/GLM paths, Cox paths, covariance-free preprocessing, marginal FDR calculations, and package-level computational helpers.

The translation is intentionally self-contained. It does not vendor or compile BLAS, LAPACK, R runtime modules, `glmnet`, `survival`, `Matrix`, or other translated R packages. A repository audit found sibling translations such as `glmnet` and `survival`, but their APIs do not implement grpreg's defining grouped/nonconvex descent and therefore are not used as computational substitutes.

Independent validation scripts under `tools/` use NumPy/SciPy only as optional test references. They are not FPM dependencies and are not required to build or use the library.

Important compatibility differences are documented in `README.md`, `API_COVERAGE.md`, and `VALIDATION.md`. RNG streams are not intended to match R exactly.
