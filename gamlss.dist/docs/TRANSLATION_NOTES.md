# Translation notes

## Scope of v0.3.0

`gamlss.dist` contains more than one hundred R files and compiled C kernels.
The Fortran project translates reusable probability and likelihood computation,
not the R family-object framework. Plotting, S3 machinery, formula/model-frame
handling, expression-valued score/EIM closures, and R diagnostics remain out of
scope.

v0.3.0 completes the major computational targets named after v0.2.0: double
binomial, ST3C, skew-normal/SST/generalized-t, ex-Gaussian/Pareto variants,
flexible gamma/NB models, and the principal remaining zero-modified
PIG/Sichel/BB/BNB/Zipf variants.

## v0.3 implementation details

- `ST3C` is exposed as the same numerical split-t law as `ST3`; upstream ST3C
  accelerates that calculation in C, which is unnecessary in native Fortran.
- `SN1` follows the Azzalini skew-normal density. Its CDF is evaluated by
  adaptive integration and its quantile by bracketed inversion.
- `SN2` uses closed-form split-normal probability calculations.
- `SST` applies the upstream standardizing moment transform to ST3.
- `GT` uses its beta-transform CDF/quantile for finite `nu`; its generalized
  error limiting case uses numerical integration/inversion.
- `exGAUS` retains the upstream normal approximation when the exponential scale
  is very small relative to the Gaussian scale.
- `DBI` computes the exact normalizing constant by finite-support log-sum-exp.
  Its fitter initializes dispersion away from exactly one because upstream's
  near-`sigma=1` binomial shortcut otherwise creates a flat numerical start.
- `PIG2` maps to PIG through the upstream transformed shape parameter.
- Zero-inflated variants retain a mixture point mass at zero; zero-adjusted
  variants assign the requested zero probability and renormalize positive mass.
- `GAF` is represented as a gamma law with
  `sigma1 = sigma * mu**(nu/2-1)`.
- `NBF` is represented as NB with local overdispersion
  `sigma1 = sigma * mu**(nu-2)`, falling back to Poisson when
  `sigma1 < 1e-4`, matching upstream.
- `ZINBF` applies the upstream zero-inflation mixture to NBF.
- Discrete PMFs return zero on noninteger observations rather than rounding.

## Fitting

`fit_gamlss` uses BFGS maximum likelihood with numerical Hessians and separate
predictors for up to four GAMLSS parameters. v0.3.0 exposes 62 supported family
constants. Fixed-denominator families use `fit_dbi`, `fit_zibb`, and `fit_zabb`.

Analytic family-specific GAMLSS score and expected-information closures are not
translated; numerical derivatives provide a portable common fitting layer.

## Provenance

Original metadata and compiled reference sources are under `upstream/`.
Selected R sources used for v0.2 and v0.3 translation checks are retained under
`upstream/reference-R-v02/` and `upstream/reference-R-v03/`; ST3 C references
used for ST3C are under `upstream/reference-src-v03/`.

The upstream package permits GPL-2 or GPL-3. This combined project is
GPL-3.0-only because several numerical support modules derive from the prior
GPL-3 VGAM Fortran translation.
