# Changelog

## 0.2.0 - 2026-07-24

- Added an exact box-index accelerated radius-neighbor search.
- Added `search_method="direct"`, `"box"`, and `"auto"` to false-neighbor, k-nearest, and Lyapunov workflows.
- Added automatic selection based on searchable sample size, embedding dimension, and scaled radius.
- Added deterministic distance-then-index tie ordering shared by both search methods.
- Added optional selected-method and exact-distance-evaluation outputs.
- Added direct-versus-box equivalence tests for neighbor indices, distances, false-neighbor calculations, and Lyapunov paths.
- Added deterministic performance-regression checks showing substantially fewer exact distance evaluations with box search.
- Added CSV application support for selecting `auto`, `box`, or `direct` search.
- Updated the API map and removed the earlier box-search omission.

## 0.1.0 - 2026-07-24

- Added modern Fortran delay embedding for univariate and multivariate series.
- Added single-scale and multi-scale correlation integrals.
- Added average mutual information.
- Added false-nearest-neighbor calculations.
- Added thresholded k-nearest-neighbor search and Kantz-style point following.
- Added Lyapunov stretching paths and linear-slope estimation.
- Added recurrence-distance and space-time separation matrices.
- Added Lorenz, Rossler, and Duffing systems.
- Added generic RK4 integration and observation callbacks.
- Added a demonstration, CSV application, and Lorenz analysis example.
- Added strict debug and optimized validation workflows.
- Preserved GPL-2.0-only licensing in `LICENSE`, `fpm.toml`, and every Fortran source file.
