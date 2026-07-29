# Changelog

## 1.0.7-fortran.1

- translated all four exported ufRisk computational functions
- implemented all six parametric model branches
- implemented normal and standardized Student-t innovations
- integrated the supplied rugarch Fortran translation
- integrated translated fracdiff and smoots numerical modules
- embedded the long-memory scale-smoothing workflow required by ufRisk
- translated the hidden ARFIMA filter-coefficient routine
- added typed options and result structures
- added explicit validation and status codes
- corrected the Christoffersen transition-count typo
- corrected double scale division in semiparametric Log-GARCH Student-t fitting
- added strict and optimized reproducible build/test script
- retained original source and dependency inputs for provenance
