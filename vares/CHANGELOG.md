# Changelog

## Unreleased

- Delegated central Student-t and F density, probability, and quantile kernels
  to `rfortran-core` while preserving the VaRES procedure names, optional
  arguments, defaults, and finite Student-t endpoint sentinels.
- Added direct log-tail and log-quantile regression checks against R values.

## 1.0.2-fortran.1

- Translated all 448 exported VaRES procedures to modern Fortran.
- Added pure elemental scalar/array APIs for 114 distribution families.
- Added self-contained special functions and inverse-CDF algorithms.
- Added transformed 96-point Gauss-Legendre expected-shortfall integration.
- Corrected documented density/CDF/quantile and probability-flag
  inconsistencies listed in `PORTING.md`.
- Added an FPM manifest, application, vector example, and five test programs.
- Added source-reference, inversion, correction, ES-reference, and core tests.
- Preserved GPL-2.0-or-later licensing, original metadata, source, and manuals.
