# Provenance

## Upstream

- R package: `mice`
- Version translated: 3.19.0
- Upstream date: 2025-12-09
- License: GPL (>= 2)
- Project: https://github.com/amices/mice
- Primary authors: Stef van Buuren and Karin Groothuis-Oudshoorn

`upstream/DESCRIPTION` and `upstream/CITATION` are copied from the supplied
`mice-master.zip`. Selected R and C++ files used as direct computational
references are retained under `upstream/R/` and `upstream/src/`.

## Translation mapping

| Fortran source | Principal upstream reference |
| --- | --- |
| `mice_native.f90` | `src/legendre.cpp` |
| `mice_matching.f90` | `src/matchindex.cpp`, `src/match.cpp`, PMM helpers |
| `mice_regression.f90` | `R/mice.impute.norm.R` and normal variants |
| `mice_impute_continuous.f90` | mean, sample, PMM, random-indicator, quadratic imputers |
| `mice_midastouch.f90` | `R/mice.impute.midastouch.R`, `R/auxiliary.R` |
| `mice_mpmm.f90` | `R/mice.impute.mpmm.R` |
| `mice_categorical.f90` | `mice.impute.logreg`, `polyreg`, augmentation helpers |
| `mice_polr.f90` | `R/mice.impute.polr.R`; proportional-odds cumulative-logit model |
| `mice_lda.f90` | `R/mice.impute.lda.R`; standard equal-covariance LDA likelihood |
| `mice_mnar.f90` | `R/mice.impute.mnar.norm.R`, `R/mice.impute.mnar.logreg.R` |
| `mice_twolevel.f90` | `2lonly.mean`, `2lonly.norm`, `2lonly.pmm` |
| `mice_2l_norm.f90` | `R/mice.impute.2l.norm.R` including `rwishart` and `symridge` |
| `mice_fcs.f90` | numerical FCS iteration in `R/sampler.R` plus univariate dispatch |
| `mice_pooling.f90` | `pool.scalar`, `barnard.rubin`, D3 numerical formulas |
| `mice_diagnostics.f90` | `md.pairs`, `md.pattern`, `flux`, `quickpred` |
| `mice_ampute.f90` | `ampute.mcar`, `ampute.continuous`, `ampute.discrete` |

## Shared dependencies

The output does not copy dependency source. It references sibling FPM packages:

```toml
rfortran-core = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

`rfortran-core` supplies the shared `dp` kind through `r_kinds`. Dense inverse
and positive-definite solve operations use `rfortran-linalg`, which in the
target repository is backed by its pinned pure-Fortran LAPACK dependency.
No system BLAS/LAPACK link directives are present in this package.

The upstream `2l.bin`/`2l.lmer` and lasso methods call external modelling
packages. This parity pass did not vendor or duplicate those engines, and it
does not declare a sibling dependency without a verified compatible callable
Fortran API.

## Translation policy

The maintained source is modern free-form Fortran. R S3 dispatch, formula and
expression parsing, data-frame/tidyverse manipulation, graphics, interactive
methods, and R object serialization are not translated. Numerical APIs accept
explicit arrays, masks, method codes, model matrices, and sensitivity-offset
matrices instead.
