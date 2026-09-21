# API coverage and compatibility

## Coverage basis

The denominator is the 48 distinct exported computational R functions from `NAMESPACE`. The following 13 exports are excluded because they are validation/configuration or presentation-only functions: `CheckData`, `CheckOptions`, `SetOptions`, `CreateBWPlot`, `CreateCovPlot`, `CreateDesignPlot`, `CreateDiagnosticsPlot`, `CreateFuncBoxPlot`, `CreateModeOfVarPlot`, `CreateOutliersPlot`, `CreatePathPlot`, `CreateScreePlot`, and `CreateStringingPlot`.

Mapped functions: **48 of 48 (100%)**. Package status remains **substantial**. The fraction records whether an exported computational function has a meaningful complete/substantial/partial mapping; it does not assert full R-interface compatibility or numerical parity across every option. The authoritative machine-readable mapping and per-function notes are in `[extra.translation]` in `fpm.toml`.

## Principal partial interfaces

- `FPCA` / `FSVD`: dense regular common-grid computational branches are translated; sparse PACE, ragged conditional-expectation scoring, and many R option branches are not.
- `MakeFPCAInputs`: the dense matrix/common-grid path is translated; R ID-vector grouping, list construction, and NA/deduplication interface branches are omitted.
- `Lwls2D` / `Lwls2DDeriv`: grid-output numerical kernels are translated; R `subset`, arbitrary `xout` point-list, and dispatch branches are omitted.
- `GetMeanCurve`, `GetCovSurface`, `GetCrCovYX`, `GetCrCovYZ`, and `GetMeanCI`: useful dense common-grid estimators are provided rather than every sparse/smoothing/options branch.
- `FPCAder`: the FPC derivative method is implemented; alternative R branches are omitted.
- `Stringing`: supported distance calculations and classical-MDS ordering are implemented; this replaces the upstream MASS `isoMDS` backend rather than copying that dependency.
- `FClust` / `kCFC`: deterministic FPCA plus k-means/reconstruction-error clustering is implemented rather than exact upstream EMCluster behavior.
- `FAM`, `MultiFAM`, and `SBFitting`: dense score-based/local-additive backfitting paths are provided without the complete upstream GAM/density/CV machinery.
- `FLM` / `FLMCI`: the scalar-response, dense-functional-predictor path and bootstrap intervals are implemented; functional-response and broader R object/formula interfaces are omitted.
- `FOptDes`, `FPCquantile`, `FCReg`, `TVAM`, and `VCAM`: computational dense/core approximations are mapped; specialized tuning, sparse/ragged, and R backend branches remain compatibility gaps.
- `WFDA`: dense normalization and deterministic monotone one-parameter warping are implemented rather than the full upstream optimization/interface surface.
- `MakeSparseGP`: finite-rank sparse GP simulation is implemented; custom R callback/covariance-function interfaces are omitted.

## Complete or near-direct numerical mappings

Direct or especially close mappings include trapezoidal integration, curve-area normalization, dynamic correlation, cross-correlation normalization, basis construction, bandwidth-neighbor calculations, support conversion, local smoothing kernels, growth-reference transformations, and several simulation helpers. Exact mapping status for each exported computational R function is recorded in `fpm.toml` and mirrored in the README table.

## RNG compatibility

`Dyn_test`, `MakeGPFunctionalData`, `MakeSparseGP`, `Sparsify`, `Wiener`, bootstrap helpers, permutation tests, and other stochastic paths use the Fortran intrinsic RNG. `seed_rng` supports deterministic Fortran tests, but no attempt is made to reproduce R's RNG stream.

## Untranslated computational exports

None. All 48 exported computational functions have a meaningful mapping, but many are intentionally marked `partial`; therefore the package is not described as a complete R replacement.
