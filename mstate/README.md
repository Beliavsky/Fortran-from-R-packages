# mstate-fortran

Modern free-format Fortran/FPM translation of the computational core of R
package `mstate` 0.3.3.

The upstream package is GPL-2-or-later. Its complete source tree is retained
under `upstream/mstate-0.3.3/` for license, attribution, and algorithm
provenance. Plotting, S3 printing, formula parsing, and R data-frame
infrastructure are not compiled.

## Implemented numerical core

- transition matrices, reachability, competing-risks and illness-death models
- `msprep`, `crprep`, landmark cutting, cross-sections and event counts
- `agmssurv.c`, Breslow/Efron hazards, covariance and Cox-backed `msfit`
- reduced-rank proportional-hazards fitting (`redrank`)
- Markov score test and wild bootstrap
- Aalen-Johansen `probtrans`, Aalen/Greenwood covariance and `ELOS`
- Nelson-Aalen, landmark AJ, grouped cumulative incidence and influence SEs
- full `mssample` forward/reset/history/censoring output modes
- relative-survival transition/state splitting and fixed covariance remapping
- translated `relsurv::expprep2` population-rate-table hazard calculations
- end-to-end `msfit.relsurv` fixed/bootstrap/both workflows
- relative-survival subject bootstrap, including original-hazard bootstrap
- bootstrap transition-probability SEs through `probtrans_bootstrap`
- `add.times` carry-forward behavior and days/years/months conversion
- HLD/HMD rate-table parsers from the integrated relsurv translation

## Relative-survival example

The native interface uses an explicit numeric rate-table covariate matrix. Its
columns must already be in the dimension order of the `ratetable_type`.

```fortran
use mstate
use relsurv_ratetable, only : ratetable_type

type(relative_msfit_type) :: rsfit

call msfit_relsurv_full(hz, ms, tr, split_trans, ratetable, xrate, rsfit, &
                        variance_mode='both', b=100, seed=123, &
                        time_format='days', info=info)
```

`rsfit%fit` contains the fixed-population Greenwood covariance mapping,
`rsfit%bootstrap` contains cumulative-hazard replicates and pointwise bootstrap
variances, and `rsfit%bootstrap_fit` stores the same point hazards with the
bootstrap variances on its diagonal covariance entries.

Bootstrap probability uncertainty is obtained with:

```fortran
call probtrans_bootstrap(rsfit%fit, rsfit%trans, rsfit%bootstrap, 0.0_dp, pt, info=info)
```

## Build

```text
fpm build
fpm test
fpm run --example demo_mstate
```

The direct validation scripts use bounds/runtime checks, strict implicit
interfaces, and floating-point traps.

## Integrated dependencies

The user-supplied survival Fortran source bundle is retained under
`vendor/survival_f90/`; the needed Cox modules are compiled directly into this
standalone package and retain LGPL-2.0-or-later SPDX declarations.

`vendor/relsurv-fortran-v0.2.0.zip` contains the complete translated relsurv
computational package. `relsurv_kinds`, `relsurv_ratetable`, and
`relsurv_parsers` are compiled here to provide the population mortality engine
and HLD/HMD rate-table loaders.

See `docs/TRANSLATION_STATUS.md` and `docs/PARITY_NOTES.md` for detailed
coverage and the remaining R-interface-only differences.
