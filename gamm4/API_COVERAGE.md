# API coverage

Coverage basis: distinct exported computational R functions. Plotting,
printing, formatting, internal helpers, and R-only dispatch/interface routines
are excluded.

Package status: **substantial**.

**1 of 1 (100%)** exported computational R functions is mapped. This percentage
measures function mapping, not complete R interface compatibility.

| R function | Fortran module | Fortran procedure(s) | Status |
| --- | --- | --- | --- |
| `gamm4` | `gamm4_fit_mod` | `gamm4_fit` | substantial |

## `gamm4`

Translated computational behavior includes:

- mgcv-style single-penalty smooth reparameterization into fixed/null and
  isotropic random components;
- ordinary grouped random effects represented by the translated lme4
  `random_term_t` API;
- Gaussian ML and REML variance profiling;
- binomial, Poisson, Gamma, inverse-Gaussian, and negative-binomial
  Laplace/PIRLS working-model fitting;
- lme4-compatible unstructured, diagonal, compound-symmetry, and AR(1)
  covariance parameterizations for ordinary random effects;
- smoothing parameters and random-effect standard deviations;
- GAM-only and full conditional linear predictors/fitted values;
- original-basis smooth coefficient reconstruction;
- upstream-style `getVb` marginal covariance calculation and EDF
  contributions.

Material differences from the R interface:

- callers provide numeric parametric designs, smooth basis/penalty objects, and
  optional lme4 random terms rather than R formulas/model frames;
- only one quadratic penalty per smooth is currently accepted; mgcv smooths
  with multiple penalties, especially `t2` terms, are rejected explicitly;
- no S3 `gam`/`merMod` object is constructed;
- no `subset`, `na.action`, factor-contrast, or environment processing is
  performed;
- the optional Python/CHOLMOD fast `getVb` route is omitted because the direct
  SPD calculation gives the same target quantity;
- variance parameters are optimized by deterministic bounded coordinate
  minimization rather than lme4's exact optimizer stack;
- exact R/lme4 warning text and convergence metadata are not reproduced.
