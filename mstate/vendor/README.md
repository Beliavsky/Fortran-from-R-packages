# Vendored dependency sources

- `survival_f90.zip`: user-supplied Fortran translation of survival.
  The Cox modules used by mstate are compiled in `src/` and retain their
  LGPL-2.0-or-later SPDX declarations.
- `relsurv-fortran-v0.2.0.zip`: complete Fortran translation of relsurv 2.3-3.
  `relsurv_kinds`, `relsurv_ratetable`, and `relsurv_parsers` are compiled in
  `src/` for mstate's relative-survival population-hazard path.
