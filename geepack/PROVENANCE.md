# Provenance

## Upstream package

- Package: `geepack`
- Upstream version translated: 1.3.13
- Upstream license: GPL (>= 3)
- CRAN package date: 2025-10-14
- Upstream authors/copyright holders: see `NOTICE.md` and `upstream/DESCRIPTION`

The exact upstream files retained in `upstream/` came from the supplied
`geepack-master.zip`; they are not edited and are excluded from the FPM source
directories.

## Translation map

| Fortran module | Principal upstream sources | Computational role |
| --- | --- | --- |
| `geepack_links` | `src/famstr.cc` | Mean links, inverse links, derivatives, variance functions |
| `geepack_correlations` | `src/famstr.cc`, `src/geesubs.cc` | Working-correlation matrices and derivatives |
| `geepack_design` | `R/genZcor.R`, `R/fixed2Zcor.R` | `genZcor`, ordinal odds design, fixed-correlation extraction |
| `geepack_gee` | `src/gee2.cc`, `src/geesubs.cc`, `src/param.cc`, `src/inter.cc` | Mean/scale/association estimating equations, sandwich covariance, influence and jackknife calculations |
| `geepack_ordinal` | `R/ordgee.R`, `src/ordgee.cc` | Clustered ordinal cumulative-response GEE and local odds ratios |
| `geepack_metrics` | `R/qic-ce.R`, `R/geese.R` | Quasi-likelihood, QIC/CIC/QICC, coefficient-comparison covariance |
| `geepack_relative_risk` | `R/relative-risk-regression.R` | COPY construction and log-binomial relative-risk GEE |
| `geepack_inference` | `R/summary.R`, `R/geeglm-anova.R` | Wald statistics and contrasts |
| `geepack_matrix` | translation support | Narrow wrappers around shared `rfortran-linalg` solves/inverses |
| `geepack_status`, `geepack` | translation support | Status codes and public facade |

The numeric link and correlation integer codes follow the upstream package.
The upstream link named `inverse` is exported here as `LINK_RECIPROCAL` because
Fortran is case-insensitive and a constant named `LINK_INVERSE` would collide
with the public procedure `link_inverse`.

## Shared dependencies

`rfortran-core` supplies `dp` from `r_kinds` and the chi-square distribution
used by Wald inference. `rfortran-linalg` supplies checked dense solves and
matrix inverses and, in the target repository, obtains LAPACK through its
pinned pure-Fortran FPM dependency. No system BLAS/LAPACK linkage is requested
by this package.

## Publications named by upstream CITATION

The retained `upstream/CITATION` should be used for full bibliographic details.
It cites work by Halekoh, Højsgaard & Yan (2006), Yan & Fine (2004), Yan
(2002), and Xu, Fine, Song & Yan (2025).
