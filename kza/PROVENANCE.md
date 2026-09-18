# Provenance

## Upstream input

Translation input: `kza-master.zip`, containing R package `kza` version 4.2.0 dated 2026-09-01.

The archive's `DESCRIPTION`, `NAMESPACE`, computational R files, and compiled C sources used to understand the numerical behavior are retained under `upstream/`.

## Files used

- `upstream/R/kza.R`: public `kz`, `kza`, `kzsv`, and `kzs` wrappers plus KZA 4.2.0 semantics.
- `upstream/R/kzft.R`: `kzft`, `kztp`, `periodogram`, and `transfer_function` computational definitions.
- `upstream/R/rlv.R`: rolling local variance implementation and boundary policies.
- `upstream/src/kz.c`: 1D/2D/3D iterated KZ moving-average kernels.
- `upstream/src/kza.c`: 1D/2D/3D KZA adaptive kernels and robust normalizer.
- `upstream/src/kzf.c`: directional-difference helper behavior.
- `upstream/src/kzsv.c`: adaptive sample-variance kernel.

## Dependency review

The target `Beliavsky/Fortran-from-R-packages` repository and its shared-module conventions were checked before implementation. This package does not require BLAS/LAPACK or another translated R package. No compatible package-specific dependency was needed for KZ/KZA filtering. A direct DFT is used for `periodogram`, so the upstream R package's FFTW system requirement is not carried into the Fortran translation.

Because no shared dependency is required, `dp` is defined once in `src/kza_kinds.f90` from `iso_fortran_env::real64`, then imported everywhere else.

## Translation approach

The Fortran code follows the corrected upstream 4.2.0 algorithms rather than older pre-4.2 behavior described in comments. In particular it preserves the full-window interpretation of `m`, iterative KZA feedback, the adaptive zero-normalizer guard, robust 99th-percentile normalization with zero fallback, fixed 2D indexing/bounds behavior, collapsed-dimension 3D averaging, and directional derivative edge handling.

Where the C sources would otherwise rely on undefined indexing for a collapsed 2D dimension, the Fortran implementation uses the intended zero directional derivative, consistent with the explicit fixes already present for the analogous 3D cases.

No upstream dependency source was vendored. The only retained upstream files are the kza package's own source and metadata for license/provenance traceability.
