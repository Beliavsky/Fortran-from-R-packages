# Changelog

## 0.4.2-fortran.1

- Translated all exported FRAPO computational routines to modern Fortran.
- Added vector and matrix overloads for return and trend operations.
- Added self-contained convex QP/LP and risk-parity solvers.
- Added long-only GMV, most-diversified, minimum-tail-dependence, and ERC portfolios.
- Added maximum-, average-, constrained-CDaR, and minimum-CDaR portfolios.
- Added typed results, status codes, FPM packaging, examples, and tests.
- Recomputed reported drawdowns from the actual running high-water mark to avoid
  nonunique LP slack-variable artifacts.
- Preserved GPL-3.0-or-later licensing, original attribution, and complete source provenance.
