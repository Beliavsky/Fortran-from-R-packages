# Rssa modern Fortran translation

This directory translates the computational core of **Rssa 1.1**, the R package for Singular Spectrum Analysis, to modern free-form Fortran organized as an FPM package. The original package sources are retained unmodified under `upstream/` for license, attribution and numerical provenance.

The implementation favors a portable dense reference path over Rssa's FFTW/C acceleration layer. It reuses the sibling `svd` translation for real/complex singular decompositions and `rfortran-linalg` for QR, least-squares, solves and general eigenvalue problems. No BLAS, LAPACK, ARPACK, FFTW or translated dependency source is vendored here.

## Build

Place this directory at the repository root next to `svd/` and `rfortran-linalg/`, then run:

```text
fpm build
fpm test
```

Convenience validation launchers are supplied for environments with FPM on `PATH`:

```text
scripts/validate.sh
```

On Windows, use `scripts\validate.bat`. The translation environment used for this release did not contain FPM; the exact direct-gfortran validation performed there is documented in `VALIDATION.md`.

## Implemented numerical families

- 1-D SSA, Toeplitz SSA, MSSA, rectangular 2-D SSA and complex SSA.
- Projection SSA with explicit row/column projector bases.
- Hankel, block-Hankel and Toeplitz dense matrix operations.
- Reconstruction, residuals, component contributions and right-vector calculation.
- Weighted norms, w-correlations and Frobenius correlations.
- Diagonal weighted-oblique SSA decomposition, I-OSSA, FOSSA, default 1-D EOSSA, and ordinary/oblique OSSA correlations.
- Linear recurrence relations and recurrent/vector/bootstrap forecasting.
- LS-ESPRIT parameter estimation for 1-D, MSSA, CSSA and rectangular 2-D SSA.
- Sequential and iterative gap filling, including 2-D/MSSA/CSSA iterative paths.
- Cadzow low-rank approximation.
- Automatic complete-link w-correlation grouping and direct-DFT spectral grouping.
- Heterogeneity matrices.

## Translation coverage

Status: **substantial**.

**122 of 122 (100.0%)** exported/registered computational functions are mapped. This fraction counts functions with meaningful numerical mappings; it does **not** mean every R option, S3 behavior, numerical acceleration path or object representation is reproduced.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `clone` | rssa_decomposition | `clone_ssa` | complete |
| `decompose` | rssa_decomposition<br>rssa_projection | `decompose_ssa`, `decompose_toeplitz`, `decompose_mssa`, `decompose_2d`, `decompose_complex`, `decompose_pssa` | substantial |
| `reconstruct` | rssa_reconstruction | `reconstruct_ssa`, `reconstruct_mssa`, `reconstruct_2d`, `reconstruct_complex` | substantial |
| `nu` | rssa_decomposition | `nu` | complete |
| `nv` | rssa_decomposition | `nv` | complete |
| `nlambda` | rssa_decomposition | `nlambda` | complete |
| `nsigma` | rssa_decomposition | `nsigma` | complete |
| `nspecial` | rssa_decomposition<br>rssa_projection | `nspecial`, `nspecial_pssa` | substantial |
| `contributions` | rssa_decomposition | `contributions` | complete |
| `calc.v` | rssa_decomposition<br>rssa_projection | `calc_v_ssa`, `calc_v_complex`, `calc_v_pssa` | substantial |
| `ssa` | rssa_decomposition<br>rssa_projection | `ssa_1d`, `ssa_mssa`, `ssa_2d`, `ssa_toeplitz`, `ssa_complex`, `decompose_pssa` | substantial |
| `wcor` | rssa_metrics | `wcor_ssa`, `wcor_default` | substantial |
| `wcor.default` | rssa_metrics | `wcor_default` | complete |
| `hmatr` | rssa_hmatr | `hmatr` | substantial |
| `wnorm` | rssa_metrics | `wnorm` | substantial |
| `new.hmat` | rssa_matrices | `new_hmat` | complete |
| `hmatmul` | rssa_matrices | `hmatmul` | complete |
| `hankel` | rssa_matrices | `hankel`, `hankel_matrix`, `hankelize_matrix` | substantial |
| `hcols` | rssa_matrices | `hcols` | complete |
| `hrows` | rssa_matrices | `hrows` | complete |
| `is.hmat` | rssa_matrices | `is_hmat` | complete |
| `new.hbhmat` | rssa_matrices | `new_hbhmat` | substantial |
| `hbhmatmul` | rssa_matrices | `hbhmatmul` | substantial |
| `hbhcols` | rssa_matrices | `hbhcols` | substantial |
| `hbhrows` | rssa_matrices | `hbhrows` | substantial |
| `is.hbhmat` | rssa_matrices | `is_hbhmat` | substantial |
| `new.tmat` | rssa_matrices | `new_tmat` | complete |
| `tmatmul` | rssa_matrices | `tmatmul` | complete |
| `tcols` | rssa_matrices | `tcols` | complete |
| `trows` | rssa_matrices | `trows` | complete |
| `is.tmat` | rssa_matrices | `is_tmat` | complete |
| `lrr` | rssa_forecast | `lrr`, `lrr_default`, `lrr_ssa`, `lrr_mssa`, `lrr_complex` | substantial |
| `roots` | rssa_forecast | `roots_lrr`, `roots_lrr_complex` | complete |
| `rforecast` | rssa_forecast | `rforecast_ssa`, `rforecast_mssa`, `rforecast_complex`, `rforecast_pssa` | substantial |
| `vforecast` | rssa_forecast | `vforecast_ssa`, `vforecast_mssa`, `vforecast_complex`, `vforecast_pssa` | substantial |
| `bforecast` | rssa_forecast | `bforecast_ssa` | partial |
| `parestimate` | rssa_parestimate | `parestimate_ssa`, `parestimate_mssa`, `parestimate_complex`, `parestimate_2d` | substantial |
| `cadzow` | rssa_cadzow | `cadzow` | substantial |
| `frobenius.cor` | rssa_metrics | `frobenius_cor` | substantial |
| `igapfill` | rssa_gapfill | `igapfill`, `igapfill_2d`, `igapfill_mssa`, `igapfill_complex` | substantial |
| `gapfill` | rssa_gapfill | `gapfill_ssa`, `gapfill_mssa_channel`, `gapfill_complex` | substantial |
| `summarize.gaps` | rssa_gapfill | `summarize_gaps`, `summarize_gaps_complex` | substantial |
| `grouping.auto` | rssa_autogroup | `grouping_auto` | substantial |
| `grouping.auto.wcor` | rssa_autogroup | `grouping_auto_wcor_ssa` | substantial |
| `grouping.auto.pgram` | rssa_autogroup | `grouping_auto_pgram_ssa` | substantial |
| `clone.ssa` | rssa_decomposition | `clone_ssa` | complete |
| `decompose.ssa` | rssa_decomposition | `decompose_ssa` | substantial |
| `decompose.toeplitz.ssa` | rssa_decomposition | `decompose_toeplitz` | substantial |
| `decompose.cssa` | rssa_decomposition | `decompose_complex` | substantial |
| `decompose.pssa` | rssa_projection | `decompose_pssa` | substantial |
| `reconstruct.ssa` | rssa_reconstruction | `reconstruct_ssa` | substantial |
| `residuals.ssa` | rssa_reconstruction | `residuals_ssa` | substantial |
| `residuals.ssa.reconstruction` | rssa_reconstruction | `residuals_ssa` | substantial |
| `calc.v.ssa` | rssa_decomposition | `calc_v_ssa` | complete |
| `calc.v.cssa` | rssa_decomposition | `calc_v_complex` | substantial |
| `calc.v.pssa` | rssa_projection | `calc_v_pssa` | substantial |
| `wcor.ssa` | rssa_metrics | `wcor_ssa` | substantial |
| `wcor.ossa` | rssa_oblique | `wcor_ossa` | substantial |
| `owcor` | rssa_oblique | `owcor_ssa` | substantial |
| `decompose.wossa` | rssa_oblique | `decompose_wossa` | substantial |
| `fossa` | rssa_oblique | `fossa_ssa` | partial |
| `fossa.ssa` | rssa_oblique | `fossa_ssa` | partial |
| `decompose.ossa` | rssa_oblique | `decompose_ossa` | complete |
| `iossa` | rssa_iterative_oblique | `iossa_ssa` | partial |
| `iossa.ssa` | rssa_iterative_oblique | `iossa_ssa` | partial |
| `eossa` | rssa_iterative_oblique | `eossa_ssa` | partial |
| `eossa.ssa` | rssa_iterative_oblique | `eossa_ssa` | partial |
| `wnorm.default` | rssa_metrics | `wnorm_default` | complete |
| `wnorm.ssa` | rssa_metrics | `wnorm_ssa` | complete |
| `wnorm.complex` | rssa_metrics | `wnorm_complex` | complete |
| `wnorm.1d.ssa` | rssa_metrics | `wnorm_ssa` | complete |
| `wnorm.cssa` | rssa_metrics | `wnorm` | complete |
| `wnorm.nd.ssa` | rssa_metrics | `wnorm_2d` | partial |
| `wnorm.toeplitz.ssa` | rssa_metrics | `wnorm_ssa` | complete |
| `wnorm.mssa` | rssa_metrics | `wnorm_mssa` | complete |
| `lrr.default` | rssa_forecast | `lrr_default` | complete |
| `lrr.1d.ssa` | rssa_forecast | `lrr_ssa` | substantial |
| `lrr.toeplitz.ssa` | rssa_forecast | `lrr_ssa` | substantial |
| `lrr.mssa` | rssa_forecast | `lrr_mssa` | substantial |
| `lrr.cssa` | rssa_forecast | `lrr_complex` | substantial |
| `forecast.1d.ssa` | rssa_forecast | `rforecast_ssa`, `vforecast_ssa`, `bforecast_ssa` | partial |
| `forecast.toeplitz.ssa` | rssa_forecast | `rforecast_ssa`, `vforecast_ssa`, `bforecast_ssa` | partial |
| `predict.1d.ssa` | rssa_forecast | `rforecast_ssa`, `vforecast_ssa`, `bforecast_ssa` | substantial |
| `predict.toeplitz.ssa` | rssa_forecast | `rforecast_ssa`, `vforecast_ssa`, `bforecast_ssa` | substantial |
| `predict.mssa` | rssa_forecast | `rforecast_mssa`, `vforecast_mssa` | partial |
| `rforecast.1d.ssa` | rssa_forecast | `rforecast_ssa` | substantial |
| `rforecast.toeplitz.ssa` | rssa_forecast | `rforecast_ssa` | substantial |
| `rforecast.mssa` | rssa_forecast | `rforecast_mssa` | partial |
| `rforecast.cssa` | rssa_forecast | `rforecast_complex` | substantial |
| `rforecast.pssa.1d.ssa` | rssa_forecast | `rforecast_pssa` | substantial |
| `vforecast.1d.ssa` | rssa_forecast | `vforecast_ssa` | substantial |
| `vforecast.toeplitz.ssa` | rssa_forecast | `vforecast_ssa` | substantial |
| `vforecast.mssa` | rssa_forecast | `vforecast_mssa` | partial |
| `vforecast.cssa` | rssa_forecast | `vforecast_complex` | substantial |
| `vforecast.pssa.1d.ssa` | rssa_forecast | `vforecast_pssa` | substantial |
| `roots.lrr` | rssa_forecast | `roots_lrr` | complete |
| `bforecast.1d.ssa` | rssa_forecast | `bforecast_ssa` | partial |
| `bforecast.toeplitz.ssa` | rssa_forecast | `bforecast_ssa` | partial |
| `parestimate.1d.ssa` | rssa_parestimate | `parestimate_ssa` | substantial |
| `parestimate.nd.ssa` | rssa_parestimate | `parestimate_2d` | partial |
| `parestimate.toeplitz.ssa` | rssa_parestimate | `parestimate_ssa` | substantial |
| `parestimate.mssa` | rssa_parestimate | `parestimate_mssa` | substantial |
| `parestimate.cssa` | rssa_parestimate | `parestimate_complex` | substantial |
| `cadzow.ssa` | rssa_cadzow | `cadzow_ssa` | substantial |
| `nspecial.ssa` | rssa_decomposition | `nspecial` | complete |
| `nspecial.pssa` | rssa_projection | `nspecial_pssa` | substantial |
| `gapfill.1d.ssa` | rssa_gapfill | `gapfill_ssa` | substantial |
| `gapfill.cssa` | rssa_gapfill | `gapfill_complex` | substantial |
| `gapfill.toeplitz.ssa` | rssa_gapfill | `gapfill_ssa` | substantial |
| `gapfill.mssa` | rssa_gapfill | `gapfill_mssa_channel` | partial |
| `igapfill.1d.ssa` | rssa_gapfill | `igapfill` | substantial |
| `igapfill.nd.ssa` | rssa_gapfill | `igapfill_2d` | partial |
| `igapfill.mssa` | rssa_gapfill | `igapfill_mssa` | substantial |
| `igapfill.cssa` | rssa_gapfill | `igapfill_complex` | substantial |
| `igapfill.toeplitz.ssa` | rssa_gapfill | `igapfill` | partial |
| `summarize.gaps.default` | rssa_gapfill | `summarize_gaps` | substantial |
| `summarize.gaps.1d.ssa` | rssa_gapfill | `summarize_gaps` | substantial |
| `summarize.gaps.cssa` | rssa_gapfill | `summarize_gaps_complex` | substantial |
| `summarize.gaps.toeplitz.ssa` | rssa_gapfill | `summarize_gaps` | substantial |
| `grouping.auto.wcor.ssa` | rssa_autogroup | `grouping_auto_wcor_ssa` | substantial |
| `grouping.auto.pgram.1d.ssa` | rssa_autogroup | `grouping_auto_pgram_ssa` | substantial |
| `grouping.auto.pgram.toeplitz.ssa` | rssa_autogroup | `grouping_auto_pgram_ssa` | substantial |

All 122 functions in the computational coverage basis have a meaningful mapping. Several mappings remain partial, so 100% mapped-function coverage does not imply complete R/S3 or option-level compatibility. See [API_COVERAGE.md](API_COVERAGE.md) for details.

## Validation

Deterministic tests cover matrix embeddings, ordinary/MSSA/2-D/projection and weighted/iterative/filtered/ESPRIT oblique decompositions, reconstruction, ordinary/oblique correlations, LRR forecasts and roots, ESPRIT, sequential and iterative gap filling, Cadzow, automatic grouping and the heterogeneity matrix. See [VALIDATION.md](VALIDATION.md) for the exact compiler/build record.

## License and provenance

Rssa declares **GPL (>= 2)** and its source headers state GPL version 2 or later. This translation uses `GPL-2.0-or-later`; the GPLv2 text is in `LICENSE`. Upstream authorship/copyright and the supplied archive checksum are recorded in `NOTICE.md` and `provenance/`.

