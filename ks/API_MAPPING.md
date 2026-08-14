# API mapping

This file maps the principal non-plotting `ks` 1.15.3 computations to the
Fortran API.  R S3 methods that only choose defaults, construct grids, print or
plot are intentionally represented by lower-level numerical procedures rather
than duplicated object wrappers.

| R `ks` functionality | Fortran API | v0.1.0 status |
|---|---|---|
| `kde`, `dkde`, `predict.kde` | `kde_model`, `fit_kde`, `kde_pdf`, `kde_logpdf` | implemented |
| `pkde`, `qkde`, `rkde` (1D) | `kde_cdf_1d`, `kde_quantile_1d`, `kde_random` | implemented |
| `kdde`, `predict.kdde` | `kdde_eval` | implemented, arbitrary derivative order |
| `kfe` | `kfe_tensor`, `kfe_1d` | implemented |
| `Hns`, `Hns.diag`, `hns` | `hns_matrix`, `hns_diag`, `hns_1d` | implemented |
| `Hlscv`/`Hucv`, `hlscv`/`hucv` | `hlscv_matrix`, `hlscv_1d`, `lscv_value` | implemented core criterion/optimizer |
| `hpi` | `hpi` | implemented staged 1D plug-in |
| `Hpi`, `Hpi.diag` | `hpi_matrix`, `hpi_diag` | implemented one-stage full-matrix plug-in; not every upstream pilot option |
| `hscv` | `hscv` | implemented 1D |
| `Hscv`, `Hbcv`, `Hnm` families | -- | deferred selector families |
| `mise.mixt`, `amise.mixt`, `ise.mixt` | `mise_normal_mixture`, `amise_normal_mixture`, `ise_normal_mixture` | implemented for density (`r=0`) mixture risk |
| `Hmise.mixt`, `Hamise.mixt` | `hmise_normal_mixture`, `hamise_normal_mixture` | implemented |
| normal-mixture density/random | `normal_mixture_pdf`, `normal_mixture_sample` | implemented |
| multivariate normal mixture modes/moments | `normal_mixture_modes`, `normal_mixture_moments` | implemented |
| Student-t mixture density/random | `student_mixture_pdf`, `student_mixture_sample` | implemented |
| normal derivative helpers | `normal_derivative`, `mvn_derivative_tensor`, `psins_1d` | implemented |
| `kda`, `predict.kda`, `compare` | `kda_model`, `fit_kda`, `predict_kda`, `classification_error`, `confusion_matrix` | implemented |
| `Hkda`/CV selector wrappers | caller-supplied `Hs` or `hns_matrix` | wrapper-specific selection deferred |
| `kcde`, `predict.kcde` | `kcde_model`, `fit_kcde`, `kcde_eval` | implemented, multivariate |
| `Hns.kcde` | `hns_matrix` | implemented equivalent Gaussian covariance rule |
| `Hpi.kcde` family | -- | deferred specialized CDF bandwidth selectors |
| `kms`, `kms.part` | `mean_shift_point`, `kms` | implemented core mean-shift/partitioning |
| `kdr` | `density_ridge_point` | pointwise ridge iteration implemented |
| `kdr.segment` | -- | deferred segmentation wrapper |
| `kcurv` | `kcurv_eval` | implemented |
| `kfs` | `kfs_eval` | implemented pointwise significance calculations |
| `histde`, `predict.histde` | `histde_1d`, `histde_2d`, `hist_predict_1d`, `hist_predict_2d` | implemented |
| `kde.boundary` beta branch | `beta_kernel2_pdf`, `boundary_kde_pdf` | implemented |
| `kcopula.de` empirical marginal branch | `pseudo_uniform_empirical`, `copula_density_empirical` | implemented |
| `kcopula` and kernel-marginal transforms | -- | deferred |
| `kdcde`/`dckde`, `reg.ucv` | `deconv_weights`, `deconv_pdf`, `reg_ucv`, `reg_ucv_value` | implemented numerical core |
| `kde.balloon` | `balloon_kde_2d` | implemented 2D |
| `kde.sp` | -- | deferred sample-point adaptive KDE |
| `ksupp` | `support_mask`, `contour_support_points`, `convex_hull_2d` | implemented density support; hull native in 2D |
| `kde.test` | `kde_two_sample_test`, `qr0` | implemented global test |
| `kde.local.test` | -- | deferred local-test object machinery |
| `kroc` | -- | deferred ROC convenience layer |
| `binning` | `linear_binning`, `grid_interpolate` | implemented generically |
| `symconv.1d`, `symconv.nd` | `symconv_1d`, `symconv_nd` | implemented |
| `vec`, `vech`, `invvec`, `invvech` | same names | implemented |
| `matrix.sqrt` | `matrix_sqrt` | implemented |
| `pre.scale`, `pre.sphere` | `pre_scale`, `pre_sphere` | implemented |
| `rowKpow` | `row_kron_power` | implemented |
| `Sdr`, `Sdrv` | `symmetrizer_matrix`, `symmetrizer_apply` | implemented |
| `Lpdiff` | `lp_grid_diff` | implemented |
| `mur`, `nur`, `nurs`, general `Qr` | -- / `qr0` for the r=0 test kernel | general high-order cumulant helpers deferred |
| plotting/color/contour display methods | -- | intentionally omitted |

## Naming note

Fortran is case-insensitive, so R names distinguished only by capitalization
cannot coexist literally.  For example, R `Hns` maps to `hns_matrix` while R
`hns` maps to `hns_1d`.
