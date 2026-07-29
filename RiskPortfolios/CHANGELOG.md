# Changelog

## 2.1.7-Fortran.1 - 2026-07-27

- Ported all executable mean, covariance, semideviation, and portfolio methods
  from RiskPortfolios 2.1.7 to modern Fortran.
- Added typed method constants and `portfolio_control`.
- Replaced R optimization dependencies with an internal constrained optimizer.
- Added BLAS/LAPACK linear algebra.
- Added independent fixed-reference numerical tests.
- Added full behavioral and constraint tests.
- Preserved GPL-2.0-or-later licensing, attribution, original metadata, and R
  source provenance.
- Documented the factor-analysis replacement and corrected source issues.
