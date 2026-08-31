# Provenance

## Input

Translation source: user-supplied archive `mitml-master.zip`.

Upstream package metadata retained in `upstream/DESCRIPTION` identifies:

- package: `mitml`
- version: `0.4-5`
- date: `2023-03-08`
- title: `Tools for Multiple Imputation in Multilevel Modeling`
- authors: Simon Grund, Alexander Robitzsch, Oliver Luedtke
- license: `GPL (>= 2)`
- `NeedsCompilation: no`

The R package contains no native C or Fortran implementation; its computational
algorithms are expressed in R and its imputation front ends delegate to the R
packages `pan` and `jomo`.

## Direct computational references

Exact upstream R files retained under `upstream/R/`:

- `internal-pool.R`
- `internal-methods-likelihood.R`
- `testEstimates.R`
- `testConstraints.R`
- `testModels.R`
- `clusterMeans.R`
- `multilevelR2.R`
- `internal-convergence.R`

The moving-average helper is defined in upstream `plot.mitml.R` and the reduced
ACF helper in `summary.mitml.R`; their small numerical formulas were translated,
while plotting/summary presentation code itself was excluded.

## Repository reuse

Before implementation, the target repository was checked for shared numerical
modules. `rfortran-core` and `rfortran-linalg` provide compatible APIs for the
shared real kind, R-compatible probability distributions, and dense linear
solves, and are therefore referenced as sibling FPM path dependencies.

The current workflow also produced separate `pan` and `jomo` translations.
They are deliberately not dependencies of this package for the licensing reason
described in `NOTICE.md`, and no source from either package is copied here.

## Translation organization

- `mitml_types.f90`: result/status data structures
- `mitml_numeric.f90`: small deterministic matrix/statistics helpers
- `mitml_pool.f90`: Rubin pooling and confidence intervals
- `mitml_tests.f90`: D1-D4 and constraint testing
- `mitml_cluster.f90`: cluster means
- `mitml_r2.f90`: multilevel R-squared/ICC formulas
- `mitml_convergence.f90`: MCMC diagnostics and autocorrelation helpers
- `mitml_likelihood.f90`: Gaussian LM/LMM log-likelihood kernels
- `mitml.f90`: public facade

No upstream R interface code is compiled by FPM.
