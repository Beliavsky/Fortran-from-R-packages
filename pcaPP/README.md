# pcaPP - modern Fortran translation

This directory translates the computational core of the R package **pcaPP
2.0-5** to modern free-form Fortran.  It is designed as a top-level package in
`Fortran-from-R-packages` and builds as a standalone FPM package without a
system BLAS/LAPACK installation.

The translation focuses on robust numerical algorithms: Qn, Kendall tau-b,
spatial L1 medians, robust projection-pursuit/grid PCA, sparse robust PCA,
covariance reconstruction, tuning scans, scaling, and the Zou simulation.
Plotting and R/S3 presentation code are intentionally omitted.

## Build

```text
fpm build
fpm test
fpm run --example robust_pca
```

No external FPM dependencies are required by the current translation.

## Main Fortran API

Use `pcapp_api`.  The principal result types are `scale_result`,
`median_result`, `pca_result`, `covariance_result`, and `tuning_result`.

For methods that are strings or function-valued arguments in R, the Fortran
API uses explicit numeric selectors:

- projected scale: `0 = sample SD`, `1 = MAD`, `2 = Qn`;
- centering in PCA/scaling paths: `-1 = none`, `0 = mean`,
  `1 = coordinate median`, `2 = spatial median`;
- `PCAproj` candidate mode: `0 = each observation`, `1 = random linear
  combinations`, `2 = random sphere directions`.

The random paths use the Fortran intrinsic RNG, so they are not intended to
match R's RNG stream bit-for-bit.

## Translation coverage

Package status: **substantial**.

**19 of 19 (100.0%)** exported computational R functions have a meaningful
Fortran mapping.  This fraction measures mapped computational functions, not
complete R compatibility.  It does not claim reproduction of R S3 classes,
plotting methods, arbitrary function-valued callbacks, optimizer bookkeeping,
all control-list fields, or exact R RNG behavior.

The presentation-only exports `PCdiagplot`, `objplot`, and `plotcov` are
excluded from the coverage denominator under the translation convention.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `PCAgrid` | `pcapp_pca` | `pca_grid` | substantial |
| `PCAproj` | `pcapp_pca` | `pca_proj` | substantial |
| `ScaleAdv` | `pcapp_pca` | `scale_adv` | substantial |
| `cor.fk` | `pcapp_api` | `cor_fk_vec_api`, `cor_fk_mat_api` | substantial |
| `covPC` | `pcapp_pca` | `cov_pc` | complete |
| `covPCAgrid` | `pcapp_pca` | `cov_pca_grid` | substantial |
| `covPCAproj` | `pcapp_pca` | `cov_pca_proj` | substantial |
| `data.Zou` | `pcapp_pca` | `data_zou` | substantial |
| `l1median` | `pcapp_l1median` | `l1median` | substantial |
| `l1median_BFGS` | `pcapp_l1median` | `l1median_bfgs` | substantial |
| `l1median_CG` | `pcapp_l1median` | `l1median_cg` | substantial |
| `l1median_HoCr` | `pcapp_l1median` | `l1median_hocr` | substantial |
| `l1median_NLM` | `pcapp_l1median` | `l1median_nlm` | substantial |
| `l1median_NM` | `pcapp_l1median` | `l1median_nm` | substantial |
| `l1median_VaZh` | `pcapp_l1median` | `l1median_vazh` | substantial |
| `opt.BIC` | `pcapp_pca` | `opt_bic` | partial |
| `opt.TPO` | `pcapp_pca` | `opt_tpo` | partial |
| `qn` | `pcapp_api` | `qn` | substantial |
| `sPCAgrid` | `pcapp_pca` | `spca_grid` | substantial |

See `docs/API_COVERAGE.md`, `docs/PROVENANCE.md`, and the function-level
entries in `fpm.toml` for compatibility notes.  Validation details are recorded
in `docs/VALIDATION.md`.
