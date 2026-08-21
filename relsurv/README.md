# relsurv-fortran

Modern free-format Fortran/FPM translation of the computational core of R package
`relsurv` 2.3-3. Plotting, R formula/model-frame plumbing, S3 printing, and other
presentation infrastructure are intentionally excluded.

Version 0.2.0 extends the initial translation with the main computational parity
targets from the upstream diagnostic, additive-model, smoothing, years-lost, and
rate-table import code.

## Build

```text
fpm build
fpm test
```

A direct GNU Fortran validation build is also possible; see `docs/VALIDATION.md`.

## Main API

```fortran
use relsurv
```

Important entry points include:

- `make_ratetable`, `pystep`, `pystep2`
- `expected_survival`, `population_survival_curve`
- `expprep2_expected`, `expprep2_summary`
- `netwei_summary`, `netfast_summary`, `population_hazard_increment`
- `rs_surv` with Pohar-Perme, Ederer I, Ederer II, and Hakulinen methods
- `cmp_rel`, `rsdiff`, `nessie_expected`
- `aalen_fit`, `aalen_fit_relative`, `aalen_fit_const`
- `rsadd_ml_rows`, `rsadd_piecewise`, `rsadd_em`, `rsadd_em_core`
- `rsadd_glm_bin`, `rsadd_glm_poisson`
- `rsadd_schoenfeld_residuals`, `rs_br`, `rs_zph`
- `rstrans_times`, `rstrans_fit`, `rsmul_fit`
- `survsplit_counting`, `inverse_time_monotone`
- `transrate`, `join_ratetables`, `transrate_hld`, `transrate_hmd`
- `epanechnikov_smooth`, `epa_smooth`, `epanechnikov_boundary_matrix`
- `years_difference`, `years_yl2013`, `years_yl2017`
- `greenwood_area_variance`, `bootstrap_column_variance`

The Cox-dependent routines use the user-supplied Fortran translation of
`survival`, embedded in this standalone source tree.

### Diagnostics

`rsadd_schoenfeld_residuals` translates the numerical residual/covariance
construction used by `residuals.rsadd`. Its output can be passed directly to
`rs_br` and `rs_zph`. Both maximum-deviation and Cramer-von Mises bridge tests,
R-style time transforms, tied-event handling, and per-event/summed covariance
scalings are available.

### Additive relative survival

`rsadd_em` implements the EM branch for unknown cause of death, including the
weighted Cox M-step, smoothed or unsmoothed excess baseline hazard, left
truncation, cause codes, automatic bandwidth comparison against the translated
Ederer-II estimator, and the missing-information covariance correction. The two
grouped GLM branches are exposed as `rsadd_glm_bin` and `rsadd_glm_poisson`.

### Rate-table file import

`transrate_hld` and `transrate_hmd` implement the computational file-parsing and
annual-probability-to-daily-hazard conversion paths for Human Life-Table and
Human Mortality Database style data. Native Fortran character arrays replace R
file/factor objects.

### Years-lost calculations

The years module includes the upstream Greenwood area variance, default years
difference, YL2013 and YL2017 integrations, curve-wise bootstrap variances, and
confidence intervals. Bootstrap replicate curves are supplied explicitly to the
Fortran API, which avoids reproducing R data-frame/report orchestration while
preserving the numerical bootstrap aggregation.

See `docs/TRANSLATION_STATUS.md` and `docs/PARITY_NOTES.md` for exact coverage and
remaining interface-level differences.

## Provenance and licensing

The complete upstream `relsurv` source is preserved under `upstream/relsurv`.
The supplied `survival_f90.zip` dependency is preserved under `vendor/`.
See `LICENSES.md` for license details.
