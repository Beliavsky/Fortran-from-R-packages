# Changelog

## 0.4.0

- Integrated `relsurv-fortran v0.2.0` rate-table machinery.
- Added `haz_function_relsurv` using translated `expprep2_summary` with left truncation.
- Added end-to-end `population_hazards_relsurv` and `msfit_relsurv`.
- Added `msfit_relsurv_full` with fixed/bootstrap/both workflows.
- Added days/years/months time conversion using the upstream 365.241-day year.
- Added `hazard_add_times` to reproduce cumulative-hazard carry-forward at requested times.
- Added `msboot_relsurv`, per-transition missing-bootstrap handling, and optional original-hazard bootstrap.
- Added `relative_bootstrap_type` and `relative_msfit_type` result containers.
- Added `probtrans_bootstrap` for bootstrap transition-probability standard errors.
- Added HLD/HMD rate-table parser modules from relsurv-fortran.
- Fixed a bounds-unsafe short-circuit assumption in the Nelson-Aalen real sorter found by bootstrap testing.
- Added dedicated rate-table, time-scaling, bootstrap, original-hazard, add-times, and probability-bootstrap tests.

## 0.3.0

- Added `crprep`, full `mssample` modes, `Cuminc` standard errors/grouping,
  multi-start-state `LMAJ` variance, path-prefix parity, and fixed-population
  relative-survival transition/covariance splitting.

## 0.2.0

- Integrated supplied survival Fortran Cox sources.
- Added stratified Cox, Cox-backed `msfit`, `redrank`, and `MarkovTest`.

## 0.1.0

- Initial modern Fortran/FPM computational translation of mstate.
