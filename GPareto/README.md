# GPareto-fortran 0.1.0

Modern Fortran/FPM translation of the computational core of the R package
**GPareto 1.1.9** (GPL-3), for Gaussian-process based multi-objective
optimization and Pareto-front uncertainty quantification.

## Main functionality

- Native Pareto dominance, non-dominated sets and exact dominated hypervolume
  for arbitrary objective dimension.
- Exact analytical bi-objective expected hypervolume improvement (EHI), ported
  from GPareto's Emmerich-Deutz C++ kernel.
- Sample-average EHI for three or more objectives and full-covariance batch
  qEHI.
- Expected maximin improvement (EMI) by sample-average approximation.
- SMS-EGO criterion.
- Stepwise uncertainty reduction (SUR) using direct Gaussian posterior
  conditioning on an integration design.
- Exact 2D/3D probability of non-domination for independent objective GPs.
- Conditional Pareto-front empirical attainment functions, Vorob'ev
  threshold/expectation/deviation (CPF).
- Conditional Pareto-set simulation and Gaussian KDE density estimation.
- Sequential `gparetoptim` and `easy_gparetoptim` workflows with optional
  covariance hyperparameter re-estimation after each evaluation.
- Bounded differential-evolution and discrete criterion optimization.
- Target-directed `get_design` optimization.
- LHS, low-discrepancy (Halton), Monte Carlo and SUR-importance integration
  designs.
- ZDT1/2/3/4/6, P1/P2, MOP2/MOP3, DTLZ1/2/3/7 and OKA1 benchmarks.

The DiceKriging numerical modules required by GPareto are vendored from the
modern Fortran DiceKriging translation under its GPL-3 license option. The
result has no R/Rcpp dependency and requires only BLAS/LAPACK.

## Build

```text
fpm build
fpm test
fpm run --example basic_gpareto
```

The manifest sets free source form, disables implicit typing and implicit
external interfaces, and links LAPACK/BLAS explicitly.

## Deliberate interface substitutions

R S3/S4/formula/data-frame machinery and all plotting/rgl code are omitted.
`rgenoud`/PSO search is replaced by bounded differential evolution. The
low-discrepancy integration-design mode uses a native Halton sequence rather
than vendoring randtoolbox's large Sobol table. Pareto-set density uses a
native diagonal Gaussian KDE instead of returning an R `ks::kde` object.

See `API_MAPPING.md` and `ALGORITHM_NOTES.md` for exact coverage and remaining
compatibility differences.
