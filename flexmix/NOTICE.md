# NOTICE and provenance

This repository is a clean Fortran translation of computational portions of the R package
`flexmix` version 2.3-21, supplied from the attached upstream source archive.

Upstream authors and contributors listed by `DESCRIPTION` are:

- Bettina Gruen — author and maintainer
- Friedrich Leisch — author
- Deepayan Sarkar — contributor
- Frederic Mortier — contributor
- Nicolas Picard — contributor

The upstream package is licensed as `GPL (>= 2)`. Many upstream R source files carry the
notice `Copyright (C) 2004-2016 Friedrich Leisch and Bettina Gruen`; those notices are
preserved verbatim in `upstream/R/`. The original `DESCRIPTION` and `NAMESPACE` are also
preserved under `upstream/`.

The Fortran files are newly written translations and carry SPDX identifier
`GPL-2.0-or-later`. They reproduce numerical ideas and formulas from the upstream package
without copying dependency source code.

## Dependency audit

The translation is intentionally self-contained and has no FPM dependencies. The shared
`Fortran-from-R-packages` repository was checked for reusable translations. In particular,
translations of `nnet` and `mclust` exist, but they are not required by the numerical API
implemented here. This avoids copying or vendoring dependency source and avoids introducing
a broad compatibility runtime solely to obtain small mixture kernels.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or translated R-package dependency source is
included or linked. Consequently there is no use of the `fortran-lapack` fork in this
package and no additional LAPACK attribution requirement.

## Scope differences

The translation accepts numeric matrices directly instead of R formulas, S4 objects,
factors, model frames, or sparse R matrix classes. The central EM algorithm, weighted and
hard classification, grouped observations, component deletion, constant and multinomial
concomitant priors, major GLM/clustering drivers including multinomial-response,
conditional-logit, robust-background, zero-inflated, fixed/shared-coefficient, penalized
elastic-net, mixed-effects, and penalized-smooth regressions, mixtures of factor analyzers,
univariate distribution mixtures, information criteria, KL divergence, relabeling,
numerical refit/score/covariance machinery, simulation, Gaussian stepwise model selection,
regression bootstrap, and bootstrap LR testing are implemented.

The newly written `flexmix_fixed.f90`, `flexmix_refit.f90`, `flexmix_bootstrap.f90`,
`flexmix_penalized.f90`, `flexmix_mixed.f90`, and `flexmix_smooth.f90` translate numerical
behavior from upstream `glmFix.R`, `refit.R`, `boot.R`, `examples.R`, `flxdist.R`,
`glmnet.R`, `lmm.R`, `lmmc.R`, `lmer.R`, and `mgcv.R` without copying dependency implementation
source. Factor-analysis and conditional-logit kernels additionally translate numerical
behavior from upstream `factanal.R` and `condlogit.R`.

`FLXMRlmmc` is translated in both its no-random (`FLXMRlmc`) and random-effect numerical
branches. The fixed-effect branch uses exact univariate truncated-normal moments; the
random-effect branch uses analytic one-dimensional truncation and a newly written
deterministic Genz/Halton quasi-Monte-Carlo kernel for multivariate truncated-normal
probabilities and moments. No `mvtnorm` source is copied or linked. R formula/model-matrix
construction (`FLXgetModelmatrix`) and the callback/list `Lapply` convenience method remain
untranslated. See `API_COVERAGE.md` and the manifest for exact mapping.
