# Translation coverage

This table describes computational parity targets in `survey` 4.5. R-only
formula/S3/graphics/database infrastructure is intentionally out of scope.

| Upstream area | Fortran status | Notes |
|---|---|---|
| `svydesign`, `svyrecvar`, multistage variance | Implemented | Numeric design object, strata/PSU/FPC, recursive stages, lonely-PSU modes. |
| `svymean`, `svytotal`, `svyvar`, `svyratio`, `svyCprod` | Implemented | Array APIs and covariance matrices. |
| `svytable`, `svycdf`, `svykappa` | Implemented | 1D/2D table primitives and design-based covariance. |
| Replicate designs / `svrVar` | Implemented core | JK1, JKn, BRR/Fay, ordinary bootstrap, and Preston rescaled multistage bootstrap (MRB). R conversion/compression conveniences omitted. |
| `calibrate`, `grake`, `rake`, `postStratify`, `trimWeights` | Implemented | Linear/GREG, raking, bounded logit, sinh. |
| `svyquantile` | Implemented core | All qrules; Taylor/Woodruff-style and replicate uncertainty. Not every historical interval option. |
| `svyglm` | Implemented core | Gaussian identity, binomial logit, Poisson log; survey sandwich covariance. |
| `svychisq` | Implemented core | Wald/adjusted-Wald and Rao-Scott paths. Weighted chi-square/F Satterthwaite and saddlepoint engines are available for model tests. |
| `svymle` | Implemented | Generic observation log-likelihood callback; `minqa` optimization, numerical Hessian, robust score sandwich when supplied. |
| `svynls` | Implemented | Generic observation-model callback with weighted Gauss-Newton and survey sandwich. |
| `svyivreg` | Implemented | Weighted 2SLS and replicate variant. |
| PPS variance | Implemented core | HT/Yates-Grundy plus Overton/Hartley-Rao/Poisson joint-inclusion helpers. |
| `svycoxph`, `svykm`, `svylogrank` | Implemented core | Uses supplied `survival` translation; robust design covariance for Cox. |
| `svysurvreg` | Implemented core | Parametric AFT via supplied `survival` translation and survey sandwich. |
| `svyttest`, `svyranktest`, `svyciprop`, contrasts | Implemented core | Array-based inference helpers. Some less-common historical CI variants remain. |
| `svycralpha`, `svyprcomp`/PCA primitives | Implemented core | Weighted correlation/covariance and PCA eigendecomposition. |
| `svyolr` | Implemented | Cumulative-link ordinal regression with logistic, probit, cloglog/Gumbel, and cauchit links; ordered cutpoints and design sandwich covariance. |
| `svyloglin` | Implemented core | Cell-probability covariance, quasi-Poisson loglinear fits, nested score/deviance comparison, weighted-eigenvalue p-values. |
| `svyfactanal` | Implemented core | ML factor analysis on survey covariance/correlation matrices, effective-n choices, Bartlett statistic, varimax. |
| two-phase / multiphase variance | Implemented numerical core | `Dcheck_strat`, multistage/subset Dcheck construction, phase Dcheck composition, phase-wise HT variance, calibration-space projection, totals and means. Formula/subset object administration omitted. |
| dual-frame / multiframe | Implemented numerical core | Two-frame constant and expected overlap weighting, totals, means, and independent-frame HT covariance. |
| score/anova model tests | Implemented numerical core | Pseudo-score, working Rao-Scott score, Wald term tests, misspecification-adjusted LRT, weighted chi-square/F Satterthwaite and saddlepoint. Formula-driven model nesting omitted. |
| `svysmooth` | External-kernel gap | Upstream `locpoly` delegates the actual smoother/bandwidth selection to `KernSmooth`; quantile smoothing delegates to `quantreg` (and splines). Those external numerical dependencies were not supplied, so the wrapper is not reproduced with a different algorithm. |
| `svysmoothArea`, `svysmoothUnit` | External dependency | Current upstream code is a thin dispatch layer to optional package `SUMMER`; no independent `survey` small-area algorithm exists here to translate. |
| CompQuadForm integration option | External dependency | Upstream exact/numerical-integration option delegates to optional `CompQuadForm`; the native Satterthwaite and saddlepoint paths are translated. |
| DBI-backed survey designs | Omitted | R/database infrastructure, not standalone numerical code. |
| graphics (`svyplot`, hist/box/QQ/coplot, panels) | Omitted | Explicitly excluded by translation request. |
| Formula/model-frame/S3 methods | Omitted | Replaced by explicit arrays/design matrices. |

The supplied MatrixExtra port is retained as reference material but is not
required by the current dense numerical implementation. Sparse-matrix parity
is a possible future performance target rather than a correctness dependency.
