# Origin and provenance

This project translates numerical routines from the R package `PerformanceAnalytics` version 2.1.0.

Original package metadata:

- Package: PerformanceAnalytics
- Version: 2.1.0
- Date: 2026-04-05
- Authors include Brian G. Peterson, Peter Carl, and package contributors listed in the original `DESCRIPTION`
- Original license: `GPL-2 | GPL-3`
- Original project: `braverock/PerformanceAnalytics`

Input archive used for the translation:

```text
PerformanceAnalytics-master.zip
SHA-256: 8893c59a72eb8373c84dfd46729d205f29750da9c297b07d9362b6cc77f71209
```

The original package contains R routines plus C kernels for co-moment estimators and helpers. This translation reimplements the self-contained numerical layers in modern Fortran and does not copy R object infrastructure or graphical code.

The original package authors did not endorse or validate this translation. Numerical differences and omissions are documented in `README.md`, `API_MAP.md`, and `VALIDATION.md`.

## Source-kernel audit for version 0.3.0

The finite-sample audit covered `R/MultivariateMoments.R` and `src/comomentsEstimators.c`. The self-contained correction surface consists of analytical VM2/VM3/VM3kstat/VM4 and target-covariance kernels; no separate lookup-table file is used for these shrinkage calculations. The Fortran translation preserves the GPL-2.0-or-later licensing and documents two corrected apparent source defects in `README.md`.
