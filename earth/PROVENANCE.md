# Provenance

## Upstream snapshot

Translation source: R package **earth 5.3.6**, from the user-supplied `earth-master.zip` snapshot.

The reference source retained in this package under `upstream/` includes:

- `DESCRIPTION` and `NAMESPACE` — package metadata and exported/S3 interfaces.
- `R/` — upstream R orchestration and method implementations used to determine computational API semantics and coverage.
- `src/` — the upstream earth-owned C/C-header numerical sources used to understand forward-pass, span, and regression behavior. The upstream `leaps.f` and `leapshdr.f` dependency-derived sources are intentionally not copied into this package.
- `tests/test.earth.R` and `tests/test.earth.Rout.save` — the upstream deterministic porting test and saved reference output used for the `trees` numerical parity regression.

Files under `upstream/` are reference/provenance material and are not part of the FPM compilation.

## Algorithmic lineage preserved from upstream

The translated numerical design follows the upstream earth package's use of multivariate adaptive regression splines (MARS), including Friedman's hinge basis construction, interaction hierarchy, minspan/endspan heuristics, forward stopping rules, effective parameter count, and GCV selection. The upstream C comments identify historical lineage through the `mda` implementation by Trevor Hastie and Rob Tibshirani and cite Jerome Friedman's MARS literature.

Upstream pruning support includes code derived through the `leaps` package and Alan Miller's subset-regression algorithms. The maintained Fortran translation currently implements backward pruning directly and does not translate or vendor `leaps.f`.

## Translation choices

- All newly maintained Fortran is free-form `.f90`.
- One real kind, `dp = real64`, is defined in `src/earth_kinds.f90` and used throughout.
- Least squares is implemented with reorthogonalized modified Gram-Schmidt, avoiding an external BLAS/LAPACK requirement.
- Forward candidates are scored by direct weighted least-squares refits instead of reproducing the upstream C incremental covariance/cache implementation.
- Candidate span placement follows earth's total-case endspan/start-span formulas and parent-active minspan progression.
- Rank-deficient hinge pairs may retain one full-rank child, paralleling the practical purpose of earth's regression-fix handling. Numerically equivalent one-child choices use a deterministic left-hinge tie break so compiler optimization does not change the pruning path.
- Backward pruning records every subset size, computes earth's effective parameter count and GCV, and selects the permitted minimum-GCV model.

See `README.md` and `API_COVERAGE.md` for intentionally omitted R-specific and advanced computational layers.
