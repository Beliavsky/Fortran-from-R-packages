# Computational API mapping

This table maps major upstream `fields` numerical workflows to the Fortran port. R plotting/S3/presentation helpers are intentionally not mapped.

| Upstream fields area | Fortran API / module | Status |
|---|---|---|
| `rdist`, `rdist.vec`, `fields.rdist.near` | `fields_rdist`, `fields_rdist_vec`, `fields_rdist_near` in `fields_distance` | translated |
| `rdist.earth`, `rdist.earth.vec`, `RdistEarth` | `fields_rdist_earth`, `fields_rdist_earth_vec` | translated |
| `compactToMat`, `addToDiagC` | `compact_to_matrix`, `add_to_diagonal` | translated |
| `Exponential`, `exp.cov`, `gauss.cov`, `double.exp`, powered exponential | scalar covariance functions and `stationary_covariance` | translated |
| `Matern`, `Matern.cor.to.range` | `matern`, `matern_cor_to_range` | translated |
| `stationary.cov`, anisotropic stationary covariance | `stationary_covariance`, `stationary_earth_covariance` | translated |
| `cubic.cov`, `RadialBasis`, `rad.cov` | `cubic_covariance`, `radial_basis`, `radial_covariance` | translated |
| `Wendland`, `stationary.taper.cov` | `wendland`, dense/sparse Wendland covariance, `stationary_taper_covariance` | translated |
| `Paciorek.cov` | `paciorek_covariance` | translated |
| `Tps.cov` | `tps_covariance`, `tps_covariance_marginal` | translated |
| `fields.mkpoly`, `makeMultiIndex`, `evlpoly*`, derivative polynomial utilities | `polynomial_basis`, `polynomial_power_table`, `make_multi_index`, evaluation/gradient routines | translated |
| `sreg`, `gcv.sreg`, `splint` | `smoothing_spline`, `smoothing_spline_gcv`, `spline_predict`, `spline_df_to_lambda` | translated |
| robust `rcss` spline computation | `robust_smoothing_spline`, `robust_spline_predict` | translated from native kernel |
| `QSreg`, `qsreg.psi`, `qsreg.sigma` | `quantile_smoothing_spline`, `qsreg_psi`, `qsreg_sigma` | translated |
| `QTps` | `quantile_thin_plate_spline` | translated |
| `Krig` | `krig_fit_covariance`, `krig_fit_stationary` | translated |
| `KrigFindLambda` / GCV / REML lambda selection | `krig_fit_stationary_gcv`, `krig_fit_stationary_reml`, profile functions | translated |
| `predict.Krig` | `krig_predict` | translated |
| `predictSE.Krig` | `krig_predict_covariance`, `krig_predict_se` | translated |
| `predictDerivative.Krig` | `krig_predict_gradient_stationary` | translated; covariance derivative evaluated numerically for general stationary families |
| `sim.Krig` | `krig_simulate` | translated |
| multiple-response Kriging | `krig_multi_fit_covariance`, `krig_multi_predict` | translated |
| `Tps` | `tps_fit`, `tps_fit_gcv`, `tps_predict` | translated |
| `mKrig` sparse solve core | `sparse_krig_fit_covariance`, `sparse_krig_fit_wendland`, `sparse_krig_fit_taper` | translated using spam |
| `predict.mKrig`, `predictSE.mKrig` | sparse prediction/covariance/SE routines | translated |
| `fastTps` | `fast_tps_fit`, `fast_tps_predict` | translated |
| `mKrigMLEGrid` | `spatial_profile_grid` | translated |
| `mKrigMLEJoint`, main numerical part of `spatialProcess` | `fit_spatial_process` | translated |
| `vgram.matrix`, `vgram` statistics | `variogram_grid` | translated, including Cressie robust form |
| `describe`, `stats.bin` | `describe_values`, `stats_bin` | translated |
| `fast.1way` | `fast_oneway` | translated |
| `make.surface.grid`, `interp.surface`, `findGridBox` | `make_surface_grid`, `bilinear_interpolate`, `find_grid_box` | translated |
| `offGridWeights1D`, `offGridWeights2D`, numerical core of `approximateCovariance2D` | `offgrid_weights_1d`, `offgrid_weights_2d` | translated |
| `multWendlandGrid` | `mult_wendland_grid` | typed wrapper around translated upstream native kernel |
| `image.smooth`, `smooth.2d` numerical convolution | `image_smooth_direct`, `image_smooth_fft` | translated |
| `interp.surface.FFT` | `fft_interp_surface` | translated for odd source-grid dimensions as upstream requires |
| `circulantEmbeddingSetup`, `circulantEmbedding` | `circulant_setup_2d`, `circulant_sample_2d` | translated, self-contained radix-2 FFT embedding |
| `sim.rf` | `simulate_random_field` | translated |
| `simSpatialData`, `sim.spatialProcess` numerical field/data generation | `simulate_spatial_data`, `conditional_field_simulation` | translated at numerical-kernel level |
| nonstationary Paciorek field simulation | `simulate_paciorek_field` | translated |
| `in.poly` | `points_in_polygon` | translated |
| `minimax.crit`, computational core of `cover.design` | `minimax_criterion`, `greedy_cover_design` | translated; greedy deterministic design replaces R's stochastic interchange orchestration |
| original `.Fortran` kernels (`css`, `rcss`, `dmaket`, `radbas`, `multrb`, `mltdrb`, `ddfind`, `inpoly`, `igpoly`, `multWendlandG`, etc.) | `fields_native` | retained and converted to free form |

## R-only / deliberately omitted surface

The following are not numerical parity targets for this Fortran library: base-R graphics, maps, color palettes, plot methods, print/summary formatting, S3 dispatch, formula/model-frame construction, package datasets, and UI-style grid/image object manipulation. Several R fast-prediction setup objects are decomposed into reusable sparse off-grid and FFT kernels rather than reproducing their R list/object layout.
