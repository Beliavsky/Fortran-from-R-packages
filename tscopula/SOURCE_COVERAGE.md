# Source coverage

## Translated numerical source groups

- `R/margins.R` -> `src/tscopula_margins.f90`
- `R/vtransforms.R` -> `src/tscopula_vtransforms.f90`
- `R/armacopula.R`, `R/sarmacopula.R`, `R/helper_vtarma.R` ->
  `src/tscopula_timeseries.f90`, `src/tscopula_compat.f90`
- `R/dvinecopula.R`, `R/dvinecopula2.R`, `R/dvinecopula3.R` ->
  `src/tscopula_paircopula.f90`, `src/tscopula_dvine.f90`
- `R/fitting_basic.R`, `R/fitting_vtscopula.R` ->
  `src/tscopula_models.f90`, `src/tscopula_compat.f90`
- `R/full_models.R`, `R/basic_objects.R` ->
  `src/tscopula_models.f90`, `src/tscopula_compat.f90`

## Omitted non-computational layers

- S4 class metadata and method registration
- plotting and graphical diagnostics
- `xts`, `zoo`, data-frame, and formula handling
- R list/name manipulation and formatted printing
- bundled example datasets

## Deliberate numerical substitutions

- External `rvinecopulib` kernels -> native Fortran pair-copula/D-vine code
- External `FKF` -> native exact ARMA Kalman filter
- External `ltsa` autocovariance -> impulse-response summation
- External `arfima` long-memory helpers -> deterministic truncation
- External `kdensity` -> Gaussian KDE empirical margin
- Joint D-vine/full-model optimization -> documented sequential/IFM fitting
